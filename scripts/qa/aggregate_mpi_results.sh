#!/usr/bin/env bash
# Summarize bench_mpi RESULT lines into a claims-ledger-ready CSV table.
#
# Usage:
#   scripts/qa/aggregate_mpi_results.sh logs/*.log > mpi-summary.csv
#   RANK_CSV=mpi-ranks.csv scripts/qa/aggregate_mpi_results.sh logs/*.log > mpi-summary.csv
#   cat run.log | scripts/qa/aggregate_mpi_results.sh

set -euo pipefail

awk '
BEGIN {
    rank_csv = ENVIRON["RANK_CSV"]
    if (rank_csv != "") {
        print "mode,lattice,rank,size,global,hostname,affinity,ompi_pml,ompi_btl,ompi_mtl,fi_provider,ucx_tls,omp_num_threads,rayon_threads,parallel" > rank_csv
    }
    print "mode,lattice,ranks,decomp,global,steps,time_s,mlups_total,mlups_per_rank,efficiency_pct,diag_calls,diag_time_s,gather_time_s,nonfinite,mass"
}
/^RESULT / {
    delete kv
    for (i = 2; i <= NF; i++) {
        split($i, p, "=")
        kv[p[1]] = p[2]
    }
    mode = kv["mode"]
    lattice = kv["lattice"]
    ranks = kv["ranks"] + 0
    key = mode "," lattice
    mlups_rank = kv["mlups_per_rank"] + 0.0
    if (!(key in base)) {
        base[key] = mlups_rank
    }
    eff = (base[key] > 0.0) ? 100.0 * mlups_rank / base[key] : 0.0
    printf "%s,%s,%d,%s,%s,%s,%s,%s,%.6f,%.2f,%s,%s,%s,%s,%s\n",
        mode, lattice, ranks, kv["decomp"], kv["global"], kv["steps"],
        kv["time_s"], kv["mlups_total"], mlups_rank, eff,
        kv["diag_calls"], kv["diag_time_s"], kv["gather_time_s"],
        kv["nonfinite"], kv["mass"]
}
/^RANK_RESULT / {
    if (rank_csv == "") {
        next
    }
    delete kv
    n = split($0, fields, " ")
    for (i = 2; i <= n; i++) {
        split(fields[i], p, "=")
        kv[p[1]] = p[2]
    }
    printf "%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n",
        kv["mode"], kv["lattice"], kv["rank"], kv["size"], kv["global"],
        kv["hostname"], csv(kv["affinity"]), kv["ompi_pml"], kv["ompi_btl"],
        kv["ompi_mtl"], kv["fi_provider"], kv["ucx_tls"], kv["omp_num_threads"],
        kv["rayon_threads"], kv["parallel"] >> rank_csv
}
function csv(s) {
    gsub(/"/, "\"\"", s)
    return "\"" s "\""
}
' "$@"
