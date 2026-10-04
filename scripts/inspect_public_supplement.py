"""Inspect a pinned author archive without executing it or digitizing figures."""

import argparse
import hashlib
import io
import json
from pathlib import Path
import zipfile

EXPECTED = "e616d985f5e71d693d57ca1ccf9a54cae3d5b43af04d4273d52cc661a84f3d5c"


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("archive")
    parser.add_argument("--out", required=True, help="new inventory JSON file")
    args = parser.parse_args(argv)
    try:
        payload = Path(args.archive).read_bytes()
        digest = hashlib.sha256(payload).hexdigest()
        if digest != EXPECTED:
            raise ValueError("author archive SHA-256 mismatch")
        embedded = []
        with zipfile.ZipFile(io.BytesIO(payload)) as archive:
            members = archive.namelist()
            for member in members:
                if member.endswith(".pptx"):
                    with zipfile.ZipFile(io.BytesIO(archive.read(member))) as presentation:
                        embedded.append({"file": member, "embedded_data_members": [name for name in presentation.namelist() if "embeddings/" in name or name.endswith((".csv", ".xlsx", ".xls", ".txt"))]})
        report = {"kind": "author-supplement-inventory-not-measurement-validation", "source_url": "https://ndownloader.figshare.com/files/15200057", "license": "CC BY 4.0", "attribution": "Sheth and Gratzl (2019), doi:10.1098/rspa.2018.0647", "archive_sha256": digest, "members": members, "embedded_data_inventory": embedded, "model_fitted": False}
        output = Path(args.out)
        output.parent.mkdir(parents=True, exist_ok=True)
        with output.open("x", encoding="utf-8") as handle:
            handle.write(json.dumps(report, indent=2) + "\n")
        print(json.dumps(report, indent=2))
    except (OSError, ValueError, zipfile.BadZipFile) as exc:
        parser.error(str(exc))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
