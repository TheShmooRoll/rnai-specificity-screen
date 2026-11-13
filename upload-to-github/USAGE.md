# rnai-specificity-screen

Design or screen *Caenorhabditis elegans* feeding RNAi inserts against the Bristol N2 WBcel235 assembly with WormBase WS298 annotation from RefSeq. The browser and the command-line pipeline use the same engine and reference.

## Run the site

From this folder, prepare the reference and serve the static site over HTTP.

```bash
python3 scripts/build_reference.py ws298-refseq.zip --download
python3 -m http.server 8000 --bind 127.0.0.1 --directory src
```

Open [http://127.0.0.1:8000](http://127.0.0.1:8000). Use a current browser with module workers and gzip decompression support. The first analysis loads about 14 MB of compressed transcript data. Genomic designs load the selected chromosome on demand. Analysis runs on your device and the app sends no insert sequence to the host.

- **Find a fragment.** Select a target gene, cDNA or genomic DNA, and all isoforms or one transcript. The initial length range is 500–1,200 bp. The range and exact-length modes both screen every possible start at the requested lengths. The results include coverage, direct flags, upstream distances when enabled, FASTA export, and a per-length comparison showing the limiting off-targets.
- **Check an insert.** Paste a single FASTA record or a raw DNA or RNA sequence. Select a target if you want the app to distinguish on-target coverage from off-targets. The `ctl-2` example loads the full annotated transcript to demonstrate the paralog flags.
- **Review a second candidate.** When available, the second fragment covers the same requested isoforms, uses disjoint reference coordinates, shares no exact 21-mers in either orientation, and shares no flagged genes or sibling isoforms with the first. The app also looks for a shorter candidate with full coverage and zero direct flags if the requested range yields no such candidate.

## Run the pipeline

Prepare the reference as above and use Node.js 20 or later. There are no npm dependencies. Command output is JSON.

```bash
npm test
node scripts/screen.mjs --help
node scripts/screen.mjs check --gene ctl-2 --insert insert.fa
node scripts/screen.mjs design --gene ctl-2 --template cdna
node scripts/screen.mjs design --gene kin-1 --template genomic --length 600
node scripts/benchmark.mjs --count 200 --seed 20261113 > benchmark.json
```

Use `--transcript` with a RefSeq accession or WormBase transcript name for one-isoform analysis. Use `--differences off` for exact 21-mers alone, `--differences 0` for exact 21-mers plus perfect 40-base windows, or `--ungapped` for a substitution-only 40-base search. `--no-upstream` skips the upstream analysis.

## Publish on GitHub Pages

The included workflow builds the pinned WS298 reference from NCBI, runs the tests, and publishes the static site. It verifies the assembly, annotation, and source checksums during the build.

1. In this repository, open **Settings → Pages**. Under **Build and deployment**, select **GitHub Actions** as the source.
2. Open **Actions → GitHub Pages** and click **Run workflow** on `main`. Once deployment succeeds, open the site link in **Settings → Pages**.

See the [GitHub Pages instructions](https://docs.github.com/en/pages/getting-started-with-github-pages/configuring-a-publishing-source-for-your-github-pages-site).

## How the results are calculated

**Direct sequence similarity.** The screen checks both insert orientations for exact 21-base matches. The optional second screen evaluates every 40-base insert window with a selectable limit of 0–3 differences. By default, each substitution, inserted base or deleted base counts as one difference, and the transcript span can be 38–42 bases. The seed indexes retain every seed position. The verifier uses banded edit distance around a surviving exact seed and can open a gap within the window. Counts are distinct matching starts in the original insert. A transcript can have many matching positions, and the JSON provides counts and representative alignments. The default is informed by the RNAi literature. See [Thakur et al.](https://doi.org/10.1534/g3.114.013052) and [Rual et al.](https://doi.org/10.1186/1471-2164-8-106).

**Candidate ranking.** Rank by requested isoform coverage, number of flagged loci or sibling isoforms, number of flagged starts, length, the minimum number of on-target 21-base starts across the requested isoforms, and annotated coding bases in the source transcript. Coverage requires at least one exact 21-base start. Compare the actual matching bases and exact starts for each isoform when judging a candidate, especially for genomic inserts with introns. For isoform-specific design, other transcripts of the same gene can appear as sibling hits. Output positions are 1-based and inclusive.

**Upstream similarity.** For each intended RNA model and off-target group, report the nearest shared window upstream of the next mapped trigger block along the mature RNA. The distance is the number of intervening bases. These distances describe sequence evidence. Efficacy and transitive silencing require experimental confirmation. The length range is an editable starting point. See [Alder et al.](https://pubmed.ncbi.nlm.nih.gov/12554873/) and [Conte et al.](https://doi.org/10.1002/0471142727.mb2603s109).

## Reference and verification

The reference contains 46,926 gene records, 56,468 searchable RNA models, and 52,089,418 RNA bases. The builder uses RefSeq RNA sequences where available and reconstructs other annotated models from exon coordinates. It omits 260 pre-miRNA records with separate mature annotations. The [reference manifest](src/data/reference.json) records the exact source and bundled asset hashes. Unknown bases cannot anchor exact matches. See [NCBI assembly GCF_000002985.6](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000002985.6/).

Rebuild with Python 3.11 or later using the NCBI Datasets archive, or have the script fetch the package and validate the pinned WS298 source hashes.

```bash
python3 scripts/build_reference.py ws298-refseq.zip
python3 scripts/build_reference.py ws298-refseq.zip --download
```

The test suite compares the seeded matcher with exhaustive alignment on synthetic fixtures, including spaced edits and a mixed indel and substitution. It also checks the `ctl-2` paralogs, the `kin-1` isoform selector and genomic model, upstream distance, and browser worker data loading.

The [200-gene benchmark](validation/benchmark-ws298.json) sampled eligible protein-coding loci reproducibly. Adding the substitution-only 40-base screen changed the top cDNA fragment in 14/200 cases. Enabling indels at the same limit changed it in 2/200 cases relative to the ungapped setting. These are changes in computational recommendations under the 500–1,200 bp search.

The Figure 4 historical Clone Mapper comparison is pending the original insert and WS235 or WS240 transcript FASTAs. Once available, run the following check. It compares matching insert starts to an independent exact scan and checks the published positive and negative pairs. See [Thakur et al., Figure 4](https://pmc.ncbi.nlm.nih.gov/articles/PMC4232539/#fig4).

```bash
node scripts/check-clone-mapper.mjs WS240 \
  sjj_B0222.9.fa F15E6.6.fa sjj_Y50E8.g.fa ZK643.8a.fa
```
