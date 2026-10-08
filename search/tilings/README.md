# Tiling certificates (2026-10-08)

Exact upper-bound certificates for 4n squares, from 2×2 starts built from the n = 241, 273, 307 records (scale ×2 and
split, then quench and hop).  Write-up: `../TILINGS.md`.

    python3 ../exact/verify_cert.py  n964.cert      # s(964)  <= 31.859493881559810546...
    python3 ../exact/verify_cert2.py n964.cert      # independent check (polygon intersection)
    sha256sum -c SHA256SUMS

Same for `n1092.cert` (s(1092) <= 33.840529687178261...) and `n1228.cert` (s(1228) <= 35.847313543154434...).
`summary.txt`: every start and variant (`../packer/tile_hop.py --summary`).
