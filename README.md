# myloasmSNPmers — SNPmer-only fork of myloasm

> [!NOTE]
> This is a fork of [myloasm](https://github.com/bluenote-1577/myloasm) stripped
> down to a single purpose: **detect SNPmers and write them out, then stop.**
>
> A SNPmer is a pair of k-mers that are identical except at their middle base,
> where both middle-base alleles pass a binomial test (minor allele is not just
> sequencing error) and a Fisher's exact strand-balance test — i.e. a
> heterozygous-SNP signature derived purely from strand-separated k-mer counts,
> without any base alignment.
>
> Running the binary performs only: read parsing → k-mer counting → SNPmer
> detection (`get_snpmers_inplace_sort`) → dump. All downstream assembly (twin
> reads, overlaps, unitig/twin graph, mapping, polishing, dereplication) has
> been removed from `main`.
>
> **Output:** `<output_dir>/snpmers.tsv` with columns:
> `split_kmer` (the k-mer with the middle base masked to 0), `mid_pos`
> (0-based middle position, `(k-1)/2`), `allele0_base`, `allele1_base`, the two
> reconstructed allele k-mers, and their strand-summed `allele0_count` /
> `allele1_count`.
>
> **Usage:** `myloasm <reads.fa|fq> -o <output_dir>` (k-mer size via `-k`,
> default 21, must be odd and < 24).
>
> **Build:** this fork drops the downstream assembly modules and their heavy
> dependencies (rust-spoa / abpoa / skani and the `src/rust-spoa`, `src/skani`
> vendored trees), so it builds from a clean checkout with a plain
> `cargo build --release` — no git submodules and no C/C++ POA/skani
> compilation required.

---

# myloasm - a new metagenome assembler for (noisy) long reads

>[!IMPORTANT]
> Documentation is hosted at https://myloasm-docs.github.io/.
>
>[Installation](https://myloasm-docs.github.io/install/), [usage](https://myloasm-docs.github.io/usage/), and more are in the documentation. 

<img src='https://raw.githubusercontent.com/myloasm-docs/myloasm-docs.github.io/refs/heads/main/docs/assets/logo-pink.svg' width='60%' />

Myloasm is a *de novo* metagenome assembler for long-read sequencing data. It takes sequencing reads and outputs polished contigs in a single command. 

See the [documentation](https://myloasm-docs.github.io/) for more information. 

### Citation

Jim Shaw, Maximillian Marin, and Heng Li. [High-resolution metagenome assembly for modern long reads with myloasm](https://www.nature.com/articles/s41587-026-03053-z). Nature Biotechnology (2026).
