# CASM-Smart-Phase

## 0.1.8

- Copy correct version of jar file to final image

## 0.1.7

- multiple pass only filters now marked as a single PASS in merge output

## 0.1.6

- Correct homozygous marked bed file reading

## 0.1.5

- No longer count non matching zygosity as within an MNV

## 0.1.4

- Fix issue in MNV merging

## 0.1.3

- generate-bed now uses column 5 (score) rather than 'name' to avoid clashes with smart-phase workings.
- smart-phase version updated to v1.2.1

## 0.1.2

- Fix bug in generate-bed output

## 0.1.1

- Add option to mark homozygous in the bed generation step.

Smart-phase doesn't phase homozygous changes, so we merge these in manually using the marked bed file at the
merged vcf step

## 0.1.0

- Add [smart-phase v1.2.0](https://github.com/paulhager/smart-phase) to Docker

## 0.0.1

- First release, no Smart-Phase, only utility python scripts
