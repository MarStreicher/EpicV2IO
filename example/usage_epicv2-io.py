import os
from pathlib import Path
from dotenv import load_dotenv

from epicv2io import BetasLoader

ROOT = Path(__file__).resolve().parent.parent
if not load_dotenv(ROOT / ".env"):
    raise FileNotFoundError(
        "Missing .env. Copy .env.example to .env and set the local paths."
    )

betas_path = Path(os.environ["BETAS_PATH"])
loader = BetasLoader(betas_path)

# probes × samples (cg only)
betas = loader.load_data()
print(betas.shape)
print(betas.iloc[:3, :3])
