// Example MongoDB schema for Tandem Repeat Domain Database

const geneExample = {
  gene_name: "TNFRSF4",
  aliases: ["CD134", "OX40"],
  chromosome: "chr1",
  location: "1q23.3",
  transcripts: [
    {
      transcript_id: "ENST00000379236",
      genomic_start: 1213735,
      genomic_end: 1214240,
      exons: [
        { exon_number: 1, start: 1213600, end: 1213734 },
        { exon_number: 2, start: 1213735, end: 1213900 },
        { exon_number: 3, start: 1213901, end: 1214240 }
      ],
      proteins: [
        {
          protein_id: "P43489",
          length: 277,
          description: "Tumor necrosis factor receptor superfamily member 4",
          status: "reviewed",
          repeats: [
            {
              repeat_type: "TNFR-Cys",
              index: 1,
              start_pos: 30,
              end_pos: 65,
              sequence: "CNCPAGHCPRHCDADRDCICGSYAGFFCNVTEGECCPDC",
              exon_mapping: ["Exon 2", "Exon 3"]
            },
            {
              repeat_type: "TNFR-Cys", 
              index: 2,
              start_pos: 66,
              end_pos: 97,
              sequence: "CRTHRSCEPGKGLVYERCCSPGTFSDQVCSPL",
              exon_mapping: ["Exon 3"]
            }
          ]
        }
      ]
    },
    {
      transcript_id: "ENST00000505014",
      genomic_start: 1213750,
      genomic_end: 1214250,
      exons: [
        { exon_number: 1, start: 1213750, end: 1213900 },
        { exon_number: 2, start: 1213901, end: 1214250 }
      ],
      proteins: [
        {
          protein_id: "A0A8V8TP56",
          length: 250,
          description: "TNF receptor superfamily member 4, isoform 2",
          status: "unreviewed",
          repeats: []
        }
      ]
    }
  ]
}
