# Biological background — GSE197576 (ITPR3 / RELB in hypoxia & colorectal cancer)

*Compiled from a deep-research literature review (fan-out web search + adversarial
fact-checking). Each claim below was verified across multiple sources.*

## Source study

This dataset is from **Moy et al., 2022, *Developmental Cell*** (PMID **35487218**):
*"Functional genetic screen identifies ITPR3/calcium/RELB axis as a driver of
colorectal cancer metastatic liver colonization."* The axis was discovered via a
genome-scale in vivo shRNA screen that identified 26 promoters of CRC liver
colonization.

## The ITPR3 → calcium → RELB axis

- **ITPR3** is a caffeine-sensitive IP3 receptor that releases calcium from the
  endoplasmic reticulum.
- This calcium release **induces expression of RELB**, a non-canonical NF-κB
  transcription factor.
- **RELB acts downstream of ITPR3 and is by itself sufficient** to drive CRC liver
  colonization.

## The hypoxia link (why this experiment exists)

> "ITPR3 and RELB drive CRC colony formation by promoting cell survival upon
> substratum detachment or **hypoxic exposure**."

So ITPR3/RELB help colorectal cancer cells **survive low oxygen** — a stress they face
during metastatic spread. That is exactly what GSE197576 probes: knock out ITPR3 or
RELB, then compare normoxia vs hypoxia.

## Therapeutic angle

- ITPR3 is **caffeine-sensitive**; pharmacologic inhibition with caffeine reduces CRC
  metastatic capacity → the axis is a **potential therapeutic target** (cell/in-vivo
  models; caffeine is a non-specific tool compound).

## Wider context

- **RELB independently acts as a CRC oncogene**: promotes growth via AKT/mTOR, confers
  5-FU chemo-resistance, signals through non-canonical NF-κB.
- **NF-κB and HIF-1α** are key, crosstalking drivers of CRC under inflammation and
  hypoxia (overlapping apoptosis / proliferation / angiogenesis / EMT pathways).
- Context-dependence caveat: the non-canonical pathway's upstream kinase **NIK can be
  tumor-attenuating**, and ITPR3 has a context-dependent tumor-suppressive arm at the
  ER–mitochondria interface — i.e. these genes are not uniformly pro-tumoral.

## How this maps onto our DESeq2 results

| Result | Interpretation |
|--------|----------------|
| ~5,000 genes change in hypoxia (Control) | the cells' broad low-oxygen survival program |
| 352 genes with ITPR3 × hypoxia interaction | candidate genes through which ITPR3 shapes hypoxic survival |
| 323 genes with RELB × hypoxia interaction | how RELB (downstream of ITPR3) shapes the hypoxic response |

## Key references

- Moy et al. 2022, Developmental Cell — https://pubmed.ncbi.nlm.nih.gov/35487218/
- RelB oncogenic role in colon cancer — https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6234565/
- NF-κB / HIF-1α crosstalk in CRC — https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5489807/
