async function getAfdbUrlAndFormat(uniprotId, preferredAccession = null) {
  // Correct AlphaFold DB endpoint (returns an array of predictions)
  const res = await fetch(`https://alphafold.ebi.ac.uk/api/prediction/${uniprotId}`);

  // If AFDB doesn’t have this UniProt ID, return nulls gracefully
  if (!res.ok) {
    console.warn(`AFDB: no prediction for ${uniprotId} (status ${res.status})`);
    return { url: null, format: null, chosenAccession: null };
  }

  const arr = await res.json();
  if (!Array.isArray(arr) || arr.length === 0) {
    console.warn(`AFDB: empty prediction list for ${uniprotId}`);
    return { url: null, format: null, chosenAccession: null };
  }

  // Choose prediction:
  // 1) If preferredAccession provided, find a prediction whose accession matches it.
  // 2) Otherwise (or if no match), fall back to the last prediction in the array.
  const pref = preferredAccession ? String(preferredAccession).trim() : "";
  let m = null;

  if (pref) {
    m = arr.find(p => {
      const acc = p.uniprotAccession || p.uniprotId || p.entryId || p.accession;
      return acc && String(acc).trim() === pref;
    }) || null;
  }

  if (!m) m = arr[arr.length - 1];

  const chosenAccession = m?.uniprotAccession || m?.uniprotId || m?.entryId || m?.accession || null;

  // Prefer mmCIF (works in 3Dmol) and only then PDB. (3Dmol does not read binaryCIF.)
  const mmcifUrl =
    m.mmcifUrl ||
    m.cifUrl ||
    m.files?.mmcif?.url ||
    m.files?.cif?.url || null;

  if (mmcifUrl) return { url: mmcifUrl, format: "mmcif", chosenAccession };

  const pdbUrl =
    m.pdbUrl ||
    m.files?.pdb?.url || null;

  if (pdbUrl) return { url: pdbUrl, format: "pdb", chosenAccession };

  console.warn(`AFDB: no mmCIF/PDB URL found for ${uniprotId}`);
  return { url: null, format: null, chosenAccession };
}

async function fetchAfdbPredictions(uniprotId) {
  const res = await fetch(`https://alphafold.ebi.ac.uk/api/prediction/${uniprotId}`);
  if (!res.ok) return [];
  const arr = await res.json();
  return Array.isArray(arr) ? arr : [];
}

function predictionAccession(p) {
  return p?.uniprotAccession || p?.uniprotId || p?.entryId || p?.accession || null;
}

export class ProteinViewer {
    constructor() {
        this.viewer = null;
        this.viewerInitialized = false;
        this.autoRotationActive = true;
        this.rotationAnimationId = null;
        this.isPageVisible = true;
        this.isViewerInViewport = true;
        this.viewerIntersectionObserver = null;
        this.handleVisibilityChangeCallback = null;
        this._isoformSelectInitialized = false;
    }

    async populateIsoformSelect(uniprotId, canonicalUniProtId = null, preferredAccession = null) {
        const sel = document.getElementById('isoformSelect');
        if (!sel) return;

        const arr = await fetchAfdbPredictions(uniprotId);

        // Build unique accession list
        const accessions = [];
        for (const p of arr) {
            const acc = predictionAccession(p);
            if (!acc) continue;
            const s = String(acc).trim();
            if (!s) continue;
            if (!accessions.includes(s)) accessions.push(s);
        }

        // If AlphaFold returns nothing, keep the existing options (don’t wipe the control)
        if (accessions.length === 0) return;

        const canon = canonicalUniProtId ? String(canonicalUniProtId).trim() : "";
        const pref = preferredAccession ? String(preferredAccession).trim() : "";

        sel.innerHTML = "";

        // IMPORTANT: only label the true canonical accession, not the currently selected one
        for (const acc of accessions) {
            const opt = document.createElement('option');
            opt.value = acc;
            opt.textContent = (canon && acc === canon) ? `${acc} (canonical)` : acc;
            sel.appendChild(opt);
        }

        // Selection priority: preferredAccession (what we intend/actually load) -> canonical -> first
        if (pref && accessions.includes(pref)) {
            sel.value = pref;
        } else if (canon && accessions.includes(canon)) {
            sel.value = canon;
        } else {
            sel.value = accessions[0];
        }

        // Attach handler once: selecting an accession forces that accession to be chosen from AFDB
        if (!this._isoformSelectInitialized) {
            sel.addEventListener('change', async () => {
                if (!this.viewer) return;

                // Clear current model
                try {
                    this.viewer.removeAllModels();
                    this.viewer.render();
                    this.viewerInitialized = false;
                } catch (_) {
                    // ignore
                }

                // Re-initialize using the selected accession as the preferredAccession,
                // while keeping the canonicalUniProtId unchanged.
                await this.initialize(uniprotId, this.processedRepeatRegions || [], canonicalUniProtId, sel.value);
            });
            this._isoformSelectInitialized = true;
        }
    }

    // canonicalUniProtId: the canonical from your data
    // preferredAccession: the user-selected accession to load (optional)
    initialize(uniprotId, repeatRegions, canonicalUniProtId = null, preferredAccession = null) {
        const proteinViewerElement = document.getElementById('proteinViewer');

        this.viewer = $3Dmol.createViewer(proteinViewerElement, {
            backgroundColor: "white",
            antialias: true,
            powerPreference: "high-performance"
        });

        return new Promise(async (resolve, reject) => {
            try {
                // Store for re-initialization from the isoform dropdown
                this.processedRepeatRegions = repeatRegions;

                // If canonical wasn't passed, try to derive it from the repeatRegions payload
                if (!canonicalUniProtId && Array.isArray(repeatRegions) && repeatRegions.length > 0) {
                    canonicalUniProtId = repeatRegions[0]?.canonicalUniProtId || null;
                }

                // Prefer user-selected accession if present; otherwise prefer canonical
                const preferred = preferredAccession || canonicalUniProtId || null;

                // Build dropdown using what we intend to load
                await this.populateIsoformSelect(uniprotId, canonicalUniProtId, preferred);

                // Ask AlphaFold which file is available, preferring the preferred accession match
                const { url, format, chosenAccession } = await getAfdbUrlAndFormat(uniprotId, preferred);

                // Sync dropdown to the actually chosen accession (API may fall back to last entry)
                const sel = document.getElementById('isoformSelect');
                if (sel && chosenAccession) {
                    const ca = String(chosenAccession).trim();
                    if (ca && Array.from(sel.options).some(o => o.value === ca)) {
                        sel.value = ca;
                    }
                }

                // Helper to try loading a given URL with a given format
                const tryLoad = (fileUrl, fmt) => {
                    $.ajax({
                        url: fileUrl,
                        success: (data) => { this.loadModel(data, repeatRegions, fmt); resolve(true); },
                        error: () => reject(new Error(`Failed to fetch ${fileUrl}`))
                    });
                };

                if (url && format) {
                    tryLoad(url, format);
                } else {
                    // Remote not available → try your local fallbacks (first mmcif, then pdb)
                    // Use chosenAccession if we have it, otherwise fall back to the input uniprotId.
                    const idForLocal = chosenAccession || canonicalUniProtId || uniprotId;
                    const localCif = `../data/AF-${idForLocal}-F1-model_v4.cif`;
                    const localPdb = `../data/AF-${idForLocal}-F1-model_v4.pdb`;

                    $.ajax({
                        url: localCif,
                        success: (data) => { this.loadModel(data, repeatRegions, "mmcif"); resolve(true); },
                        error: () => {
                            $.ajax({
                                url: localPdb,
                                success: (data) => { this.loadModel(data, repeatRegions, "pdb"); resolve(true); },
                                error: () => {
                                    console.log(`No 3D model available for ${idForLocal}`);
                                    this.showLoadingError(`No 3D structure available for ${idForLocal}`);
                                    reject(new Error(`No 3D model available for ${idForLocal}`));
                                }
                            });
                        }
                    });
                }
            } catch (err) {
                console.error(err);
                this.showLoadingError("Error contacting AlphaFold");
                reject(err);
            }
        });
    }

    loadModel(data, repeatRegions, format = "pdb") {
        try {
            this.hideLoadingIndicator();

            // Use the correct parser: "mmcif" or "pdb"
            const model = this.viewer.addModel(data, format);

            this.highlightRepeats(repeatRegions);

            if (repeatRegions.length > 0) {
                const firstRepeat = repeatRegions[0].start;
                const lastRepeat = repeatRegions[repeatRegions.length - 1].end;
                this.viewer.zoomTo({ resi: `${firstRepeat}-${lastRepeat}` });
            }

            this.viewer.rotate(30, { x: 1 });
            this.viewer.rotate(20, { y: 1 });

            if (typeof this.viewer.enableSlabbing === 'function') {
                this.viewer.enableSlabbing();
            }

            this.viewer.render();
            this.viewerInitialized = true;

            setTimeout(() => {
                this.highlightRepeats(repeatRegions);
                this.startAutoRotation();
                this.setupManualRotationDetection();
            }, 100);
        } catch (error) {
            console.error("Error in loadModel function:", error);
            this.showLoadingError("Error loading protein structure");
        }
    }

    highlightRepeats(repeatRegions) {
        if (!this.viewer || !this.viewerInitialized) return;
        
        try {
            this.viewer.removeAllSurfaces();
            this.viewer.removeAllShapes();
            this.viewer.removeAllLabels();
            
            this.viewer.setStyle({}, {
                cartoon: {color: '0xcccccc', opacity: 0.5},
                surface: {opacity: 0.6, color: '0xdddddd'}
            });
            
            repeatRegions.forEach(repeat => {
                const selection = {resi: `${repeat.start}-${repeat.end}`};
                
                this.viewer.setStyle(selection, {
                    cartoon: {color: repeat.colorHex, opacity: 1.0},
                    surface: {opacity: 0.8, color: repeat.colorHex}
                });
                
                this.viewer.addLabel(repeat.label, {
                    position: {resi: Math.floor((repeat.start + repeat.end) / 2)},
                    backgroundColor: 'black',
                    backgroundOpacity: 0.7,
                    fontColor: 'white',
                    fontSize: 12
                });
            });
            
            this.viewer.zoomTo();
            this.viewer.render();
        } catch (error) {
            console.error("Error in highlightRepeats:", error);
        }
    }

    showStandardView(repeatRegions) {
        if (!this.viewer || !this.viewerInitialized) return;
        
        try {
            this.viewer.removeAllSurfaces();
            this.viewer.removeAllShapes();
            this.viewer.removeAllLabels();
            this.viewer.setStyle({}, {});
            this.viewer.setStyle({}, {cartoon: {colorscheme: 'chainHetatm', thickness: 0.8}});
            
            repeatRegions.forEach(repeat => {
                this.viewer.addStyle({resi: `${repeat.start}-${repeat.end}`}, {
                    cartoon: {opacity: 1.0, thickness: 1.2}
                });
            });
            
            this.viewer.zoomTo();
            this.viewer.render();
        } catch (error) {
            console.error("Error in showStandardView:", error);
        }
    }

    zoomToRepeats(repeatRegions) {
        if (!this.viewer || !this.viewerInitialized || !repeatRegions.length) return;
        
        try {
            const firstRepeat = repeatRegions[0].start;
            const lastRepeat = repeatRegions[repeatRegions.length - 1].end;
            
            this.viewer.zoomTo({resi: `${firstRepeat}-${lastRepeat}`});
            this.viewer.rotate(30, {x: 1});
            this.viewer.rotate(20, {y: 1});
            this.viewer.render();
        } catch (error) {
            console.error("Error in zoomToRepeats:", error);
        }
    }

    showExonBoundaries(exonData) {
        if (!this.viewer || !exonData?.exons?.length) return;
        
        this.viewer.removeAllShapes();
        this.viewer.removeAllLabels();
        
        const positionMap = new Map();
        
        exonData.exons.forEach(exon => {
            if (!exon.repeat_start || !exon.repeat_end) return;
            
            if (!positionMap.has(exon.repeat_start)) {
                positionMap.set(exon.repeat_start, { labels: [] });
            }
            positionMap.get(exon.repeat_start).labels.push(`E${exon.exon_number} start`);
            
            if (!positionMap.has(exon.repeat_end)) {
                positionMap.set(exon.repeat_end, { labels: [] });
            }
            positionMap.get(exon.repeat_end).labels.push(`E${exon.exon_number} end`);
        });
        
        for (const [position, info] of positionMap.entries()) {
            const atoms = this.viewer.getModel().selectedAtoms({resi: position, atom: "CA"});
            if (atoms.length === 0) continue;
            
            const atom = atoms[0];
            let sphereColor = '0x00FF00';
            
            const hasStart = info.labels.some(l => l.includes('start'));
            const hasEnd = info.labels.some(l => l.includes('end'));
            
            if (hasStart && hasEnd) {
                sphereColor = '0xFF00FF';
            } else if (hasEnd) {
                sphereColor = '0xFF0000';
            }
            
            this.viewer.addSphere({
                center: {x: atom.x, y: atom.y, z: atom.z},
                radius: 1.5,
                color: sphereColor,
                wireframe: true,
                linewidth: 3,
                alpha: 0.9
            });
            
            info.labels.forEach((label, idx) => {
                const offset = (idx - (info.labels.length - 1) / 2) * 2;
                const bgColor = label.includes('start') ? '0x004400' : '0x440000';
                
                this.viewer.addLabel(label, {
                    position: {x: atom.x, y: atom.y + offset, z: atom.z},
                    backgroundColor: bgColor,
                    backgroundOpacity: 0.8,
                    fontColor: 'white',
                    fontSize: 12
                });
            });
        }
        
        this.viewer.render();
    }

    hideExonBoundaries() {
        if (!this.viewer) return;
        this.viewer.removeAllShapes();
        this.viewer.removeAllLabels();
        this.viewer.render();
    }

    startAutoRotation() {
        if (!this.viewer || !this.viewerInitialized) return;
        
        this.autoRotationActive = true;
        let lastFrameTime = performance.now();
        const ROTATION_SPEED = 5;
        let accumulatedRotation = 0;
        
        this.setupIntersectionObserver();
        this.setupVisibilityListener();
        
        const rotateStep = (timestamp) => {
            if (!this.autoRotationActive || !this.viewer || !this.isPageVisible || !this.isViewerInViewport) {
                this.rotationAnimationId = null;
                return;
            }
            
            const deltaTime = Math.min((timestamp - lastFrameTime) / 1000, 0.1);
            lastFrameTime = timestamp;
            
            const rotationAngle = ROTATION_SPEED * deltaTime;
            accumulatedRotation += rotationAngle;
            
            if (accumulatedRotation >= 0.5) {
                this.viewer.rotate(accumulatedRotation, {y: 1});
                this.viewer.render();
                accumulatedRotation = 0;
            }
            
            this.rotationAnimationId = requestAnimationFrame(rotateStep);
        };
        
        if (this.isPageVisible && this.isViewerInViewport && !this.rotationAnimationId) {
            this.rotationAnimationId = requestAnimationFrame(rotateStep);
        }
    }

    stopAutoRotation() {
        this.autoRotationActive = false;
        
        if (this.rotationAnimationId) {
            cancelAnimationFrame(this.rotationAnimationId);
            this.rotationAnimationId = null;
        }
        
        if (this.handleVisibilityChangeCallback) {
            document.removeEventListener('visibilitychange', this.handleVisibilityChangeCallback);
            this.handleVisibilityChangeCallback = null;
        }
        
        if (this.viewerIntersectionObserver) {
            this.viewerIntersectionObserver.disconnect();
            this.viewerIntersectionObserver = null;
        }
    }

    setupIntersectionObserver() {
        if ('IntersectionObserver' in window) {
            const viewerElement = document.getElementById('proteinViewer');
            if (this.viewerIntersectionObserver) {
                this.viewerIntersectionObserver.disconnect();
            }
            this.viewerIntersectionObserver = new IntersectionObserver((entries) => {
                this.isViewerInViewport = entries[0].isIntersecting;
            }, { threshold: 0.1 });
            this.viewerIntersectionObserver.observe(viewerElement);
        }
    }

    setupVisibilityListener() {
        if (this.handleVisibilityChangeCallback) {
            document.removeEventListener('visibilitychange', this.handleVisibilityChangeCallback);
        }
        this.handleVisibilityChangeCallback = () => {
            this.isPageVisible = document.visibilityState === 'visible';
        };
        document.addEventListener('visibilitychange', this.handleVisibilityChangeCallback);
    }

    setupManualRotationDetection() {
        const viewerElement = document.getElementById('proteinViewer');
        viewerElement.addEventListener('mousedown', () => this.stopAutoRotation());
        viewerElement.addEventListener('touchstart', () => this.stopAutoRotation());
    }

    hideLoadingIndicator() {
        const loadingIndicator = document.querySelector('.loading-indicator');
        if (loadingIndicator) {
            loadingIndicator.style.display = 'none';
        }
    }

    showLoadingError(message) {
        const loadingIndicator = document.querySelector('.loading-indicator');
        if (loadingIndicator) {
            loadingIndicator.textContent = message;
        }
    }
}
