"""Fetch bounded public source-page snapshots; never label them measurement data."""

import argparse
import hashlib
import json
from pathlib import Path
import urllib.error
import urllib.request


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", required=True, help="new local snapshot directory")
    args = parser.parse_args(argv)
    output = Path(args.out)
    if output.exists():
        parser.error("output already exists")
    catalog_path = Path(__file__).resolve().parents[1] / "docs/public-source-qualification.json"
    catalog_bytes = catalog_path.read_bytes()
    catalog = json.loads(catalog_bytes)
    output.mkdir(parents=True)
    (output / "qualification_catalog_at_probe.json").write_bytes(catalog_bytes)
    probes = []
    for source in catalog["sources"]:
        probe = {"id": source["id"], "url": source["url"], "measurement_data_acquired": False}
        try:
            request = urllib.request.Request(source["url"], headers={"User-Agent": "offline-research-source-audit/1.0"})
            with urllib.request.urlopen(request, timeout=30) as response:
                payload = response.read(2_000_001)
                if len(payload) > 2_000_000:
                    raise ValueError("page exceeds 2 MB audit limit")
                probe.update(http_status=response.status, final_url=response.url, content_type=response.headers.get("Content-Type"), snapshot_sha256=hashlib.sha256(payload).hexdigest(), snapshot_bytes=len(payload))
            snapshot = source["id"] + ".snapshot"
            (output / snapshot).write_bytes(payload)
            probe["snapshot_file"] = snapshot
        except (OSError, ValueError, urllib.error.HTTPError) as exc:
            probe["fetch_error"] = str(exc)
        probes.append(probe)
    report = {"schema_version": 1, "kind": "source-page-access-probe-not-measurement-data", "catalog_sha256": hashlib.sha256(catalog_bytes).hexdigest(), "probes": probes,
              "limits": ["Page hashes identify these snapshots only; dynamic HTML can change on a repeat fetch.", "HTTP success can include a robot-check page and does not establish measurement-file availability.", "The curated catalog records source eligibility; this script probes page access and hashes bytes only."]}
    (output / "access_probe.json").write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(report, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
