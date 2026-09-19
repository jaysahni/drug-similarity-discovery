import base64
import hashlib
import io
import json
import ssl
import tarfile
import urllib.request
from pathlib import Path

URL = "https://registry.npmjs.org/3dmol/-/3dmol-2.4.2.tgz"
INTEGRITY = "7R6uFaF5++s9EpO+R0a3bs6igi5idI+sgrVhtozaZlOGHg6oourYkAlAwYxWRaeBGlNcAd+oEIhbhJYhNxGyZA=="


def main():
    ca = Path("/etc/ssl/cert.pem")
    context = ssl.create_default_context(cafile=str(ca) if ca.exists() else None)
    with urllib.request.urlopen(URL, timeout=60, context=context) as response:
        archive = response.read(20_000_000)
    if base64.b64encode(hashlib.sha512(archive).digest()).decode() != INTEGRITY:
        raise ValueError("Pinned 3Dmol archive integrity mismatch")
    with tarfile.open(fileobj=io.BytesIO(archive), mode="r:gz") as package:
        javascript = package.extractfile("package/build/3Dmol-min.js").read()
        licenses = [m for m in package.getmembers() if m.isfile() and m.name.lower() in ("package/license", "package/license.txt", "package/license.md")]
        if len(licenses) != 1:
            raise ValueError("Expected exactly one upstream license")
        license_text = package.extractfile(licenses[0]).read()
    out = Path(__file__).resolve().parent / "vendor"
    out.mkdir(exist_ok=True)
    (out / "3Dmol-min.js").write_bytes(javascript)
    (out / "3Dmol-LICENSE.txt").write_bytes(license_text)
    (out / "provenance.json").write_text(json.dumps({"package": "3dmol", "version": "2.4.2", "source": URL, "archive_sha512": INTEGRITY, "javascript_sha256": hashlib.sha256(javascript).hexdigest()}, indent=2))
    print("Vendored 3Dmol 2.4.2 with verified npm integrity and upstream license")


if __name__ == "__main__":
    main()
