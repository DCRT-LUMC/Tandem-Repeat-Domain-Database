from fastapi import FastAPI

app = FastAPI(title="TReXomeDB API")


@app.get("/")
def root():
    return {"message": "TReXomeDB API is running"}
