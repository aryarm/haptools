# Changelog

## [0.5.0](https://github.com/CAST-genomics/haptools/compare/v0.4.2...v0.5.0) (2024-12-11)


### Features

* a new `Haplotypes.merge()` method ([#256](https://github.com/CAST-genomics/haptools/issues/256)) ([cf4ccb2](https://github.com/CAST-genomics/haptools/commit/cf4ccb2822d6b1fd42ca1208ad8c78a32fa107e1))
* allow for specifying the order of population labels in `Breakpoints.encode()` ([#262](https://github.com/CAST-genomics/haptools/issues/262)) ([09a916a](https://github.com/CAST-genomics/haptools/commit/09a916a3bf4c7e3494507af1d13673e2b78c8ff3))
* allow multiple variants in `pearson_corr_ld()` ([#258](https://github.com/CAST-genomics/haptools/issues/258)) ([808c31e](https://github.com/CAST-genomics/haptools/commit/808c31ed12268089ee960c0f04d5cd1c916d718c))
* support for Apple M1 silicon ([#255](https://github.com/CAST-genomics/haptools/issues/255)) ([3f055ac](https://github.com/CAST-genomics/haptools/commit/3f055acdb6f4effc58fde6dcff993654dbf73985))
* support for numpy 2.0 ([#250](https://github.com/CAST-genomics/haptools/issues/250)) ([69192d0](https://github.com/CAST-genomics/haptools/commit/69192d0a5fc7dcd5ca15260b4505306ff742bf1b))


### Bug Fixes

* 1) improve error in `transform` when providing an empty hap file and a `--region` and 2) allow for calling `write()` on Genotypes objects without variants ([#264](https://github.com/CAST-genomics/haptools/issues/264)) ([4e84178](https://github.com/CAST-genomics/haptools/commit/4e8417899d01b2e14bcb4c96971337817637caa0))
* allow `simphenotype` to accept TR PGENs without `--repeats` ([#263](https://github.com/CAST-genomics/haptools/issues/263)) ([16a84d1](https://github.com/CAST-genomics/haptools/commit/16a84d1fb5e1dee335b13fa0628e9abf4e53453a))
* issue warnings instead of errors when trying to check the header of hap files and issue error when output of transform is not provided to simphenotype ([#254](https://github.com/CAST-genomics/haptools/issues/254)) ([0226653](https://github.com/CAST-genomics/haptools/commit/02266539cb5d242497658e9fffec3aa7139fe8ee))
* load an empty array if there are no variants in `GenotypesVCF.read` ([#257](https://github.com/CAST-genomics/haptools/issues/257)) ([9eedc3a](https://github.com/CAST-genomics/haptools/commit/9eedc3afca202c0dd83bcce9cef34f339df83118))


### Documentation

* creating an `.hanc` file from a `.bp` file ([#260](https://github.com/CAST-genomics/haptools/issues/260)) ([b6dcd0c](https://github.com/CAST-genomics/haptools/commit/b6dcd0c44a524e2bdae2bad8f6e31fe25dee752e))

## [0.4.2](https://github.com/CAST-genomics/haptools/compare/v0.4.1...v0.4.2) (2024-06-14)


### Bug Fixes

* mention that versions must be synced in CI/CD ([35e14c0](https://github.com/CAST-genomics/haptools/commit/35e14c0da07325985c72c0919c02d245f0a3c635))

## [0.4.1](https://github.com/CAST-genomics/haptools/compare/v0.4.0...v0.4.1) (2024-06-13)


### Bug Fixes

* disable numpy 2.0 since our dependencies don't yet work with it ([#246](https://github.com/CAST-genomics/haptools/issues/246)) ([7f646b2](https://github.com/CAST-genomics/haptools/commit/7f646b2e193cbfe27b8bed07f3d10976a53b0d50))


### Documentation

* move google drive images to github ([#240](https://github.com/CAST-genomics/haptools/issues/240)) ([e8b60c0](https://github.com/CAST-genomics/haptools/commit/e8b60c0a759a4fc4f03b610fb74dbae69eca2e2d))

## [0.4.0](https://github.com/CAST-genomics/haptools/compare/v0.3.0...v0.4.0) (2024-01-14)


### Features

* a new `GenotypesPLINKTR` class for reading TRs from PGEN files ([#222](https://github.com/CAST-genomics/haptools/issues/222)) ([3c7abe6](https://github.com/CAST-genomics/haptools/commit/3c7abe6922e8d6953debf3b0d6a02dee9610bbda))
* allow multiallelic variants in `transform` ([#232](https://github.com/CAST-genomics/haptools/issues/232)) ([371415c](https://github.com/CAST-genomics/haptools/commit/371415cb9007522221ee800ab0c9570703aae6d4))
* support for python 3.11 ([#207](https://github.com/CAST-genomics/haptools/issues/207)) ([8e01ed4](https://github.com/CAST-genomics/haptools/commit/8e01ed449e47eda6cf032504c9536b053f97d588))


### Bug Fixes

* `UnboundLocalError` arising from headerless `.hap` files ([#229](https://github.com/CAST-genomics/haptools/issues/229)) ([a499b0c](https://github.com/CAST-genomics/haptools/commit/a499b0cd294c7709f866c330f7ce7431320acfbd))
* bug where `Phenotypes.subset(inplace=True)` would raise an AttributeError ([#226](https://github.com/CAST-genomics/haptools/issues/226)) ([cff6d9b](https://github.com/CAST-genomics/haptools/commit/cff6d9b0082a6170c7b0873b25445fd0915f9aab))
* convert `samples` argument in `Genotypes.read` into a set and fix `tr_harmonizer` bug arising when TRTools is also installed ([#225](https://github.com/CAST-genomics/haptools/issues/225)) ([06cc273](https://github.com/CAST-genomics/haptools/commit/06cc273fb063d73ec65071c3807c76bf63f0d448))
* Not having 23 chromosomes in genotype blocks when 23 chromosomes listed in centromere file resulted in Value Error ([#234](https://github.com/CAST-genomics/haptools/issues/234)) ([ef36798](https://github.com/CAST-genomics/haptools/commit/ef3679838e0e0193eefad2cacc37fc4d1b26715d))


### Documentation

* fix example of `.blocks.det` to `.hap` conversion in API docs ([#236](https://github.com/CAST-genomics/haptools/issues/236)) ([1ed9139](https://github.com/CAST-genomics/haptools/commit/1ed9139b5425c754256e1fc981259efbbde35d62))
* handle whitespace in blocks2hap example ([#237](https://github.com/CAST-genomics/haptools/issues/237)) ([bbdacf8](https://github.com/CAST-genomics/haptools/commit/bbdacf8f53522237ea846ca1a8f2c59c978a1d7d))

## [0.3.0](https://github.com/CAST-genomics/haptools/compare/v0.2.1...v0.3.0) (2023-06-02)


### Features

* `.snplist` input to `simphenotype` ([#217](https://github.com/CAST-genomics/haptools/issues/217)) ([cb18970](https://github.com/CAST-genomics/haptools/commit/cb18970ca65ff8ee85dc83ea68c63af3c65fda59))
* Added ability to read tandem repeats with GenotypesTR  ([#204](https://github.com/CAST-genomics/haptools/issues/204)) ([6257264](https://github.com/CAST-genomics/haptools/commit/6257264c655752cac1324ebdd0387f1207c433d6))
* Added ability to set vcftype for reading str files ([#214](https://github.com/CAST-genomics/haptools/issues/214)) ([0d734cd](https://github.com/CAST-genomics/haptools/commit/0d734cdd9235dca128af6d53fbf043f0f17511ee))
* Clump ([#211](https://github.com/CAST-genomics/haptools/issues/211)) ([3740ec1](https://github.com/CAST-genomics/haptools/commit/3740ec108075594d3b43c81840600ac59283a940))
* do not require sorting `.hap` lines by line type ([#208](https://github.com/CAST-genomics/haptools/issues/208)) ([f221397](https://github.com/CAST-genomics/haptools/commit/f2213979f66b954fa21e4ed8b953426f1736d9f8))
* new `Phenotypes.subset()` method ([#203](https://github.com/CAST-genomics/haptools/issues/203)) ([c5594d9](https://github.com/CAST-genomics/haptools/commit/c5594d9b5d6455bd9a231f218610153117b2eec7))
* simphenotype `--environment` option ([#216](https://github.com/CAST-genomics/haptools/issues/216)) ([bf69147](https://github.com/CAST-genomics/haptools/commit/bf69147f1e149d07d52adb6d407e5ec8a5e91ed3))
* Simphenotype and Index Repeat Support ([#209](https://github.com/CAST-genomics/haptools/issues/209)) ([9e2ffe1](https://github.com/CAST-genomics/haptools/commit/9e2ffe1f432459ae424bf3c2cb5cdf4f78fb98dc))


### Bug Fixes

* `Covariates.__init__` after updates to parent class ([#206](https://github.com/CAST-genomics/haptools/issues/206)) ([ce2337b](https://github.com/CAST-genomics/haptools/commit/ce2337bfb5f295942dfee2e6d1bd482a440c1d5e))
* Added logic to finding all coord files ([#201](https://github.com/CAST-genomics/haptools/issues/201)) ([be1d992](https://github.com/CAST-genomics/haptools/commit/be1d992d9a3d7651b05f66c968df354c2f71747e))
* check missing to check for 254 ([#213](https://github.com/CAST-genomics/haptools/issues/213)) ([afeab85](https://github.com/CAST-genomics/haptools/commit/afeab85bb9b70e95ca937c253a57e5fa98599b22))
* explicitly ignore repeats in the `ld` command ([#218](https://github.com/CAST-genomics/haptools/issues/218)) ([b9d0da1](https://github.com/CAST-genomics/haptools/commit/b9d0da13b7ff089ca31b8101a3b79de6806db224))
* GenotypesTR to properly load repeat count instead of GT ([#212](https://github.com/CAST-genomics/haptools/issues/212)) ([93a4eb2](https://github.com/CAST-genomics/haptools/commit/93a4eb285dc2d556cc47033750723b8313a17934))

## [0.2.1](https://github.com/CAST-genomics/haptools/compare/v0.2.0...v0.2.1) (2023-03-22)


### Bug Fixes

* NoneType error in `Haplotypes.__iter__` ([#197](https://github.com/CAST-genomics/haptools/issues/197)) ([aade751](https://github.com/CAST-genomics/haptools/commit/aade751001fb7f008382d246663dfee68886d6c6))
* precision of phenotypes written to pheno file ([#199](https://github.com/CAST-genomics/haptools/issues/199)) ([a397c96](https://github.com/CAST-genomics/haptools/commit/a397c964d80c6a1a947458a1f1b0393974ae7102))

## [0.2.0](https://github.com/CAST-genomics/haptools/compare/v0.1.3...v0.2.0) (2023-03-06)


### Features

* `Phenotypes.check missing()` method ([#191](https://github.com/CAST-genomics/haptools/issues/191)) ([621fc62](https://github.com/CAST-genomics/haptools/commit/621fc624d8b58223349767139ddd93bab6d622d3))
* Sampling without replacement option for simgenotype ([#194](https://github.com/CAST-genomics/haptools/issues/194)) ([85bd494](https://github.com/CAST-genomics/haptools/commit/85bd4946fe70e72d8e18175cc9b19acbd9e39299))


### Bug Fixes

* Fixed faulty coord file parsing logic ([#196](https://github.com/CAST-genomics/haptools/issues/196)) ([f9819b1](https://github.com/CAST-genomics/haptools/commit/f9819b15416eb1c596f8e84167f1336f6885fdd6))
* regression in multiallelic support for `simgenotype` ([#195](https://github.com/CAST-genomics/haptools/issues/195)) ([b57f91f](https://github.com/CAST-genomics/haptools/commit/b57f91fcd519845de23ce731648278bec6e3d056))


### Documentation

* describe how to add to our CLI ([#193](https://github.com/CAST-genomics/haptools/issues/193)) ([2056c7e](https://github.com/CAST-genomics/haptools/commit/2056c7e8488fdb5b2c3e47c10b98e8428fa0f07a))

## [0.1.3](https://github.com/CAST-genomics/haptools/compare/v0.1.2...v0.1.3) (2023-02-10)


### Bug Fixes

* build workflow for comment bot ([#186](https://github.com/CAST-genomics/haptools/issues/186)) ([ebcf6cd](https://github.com/CAST-genomics/haptools/commit/ebcf6cd508701951ce2f1785eb3bcfee6aa69041))
* comment bot in build pipeline ([#184](https://github.com/CAST-genomics/haptools/issues/184)) ([b3fb839](https://github.com/CAST-genomics/haptools/commit/b3fb83939c6ec9141dc9c0d9b89836cfe0fee394))


### Documentation

* warn about pysam pip issue ([#187](https://github.com/CAST-genomics/haptools/issues/187)) ([bda0157](https://github.com/CAST-genomics/haptools/commit/bda015751dc5e5aa5d3672e77737d837cd116782))

## [0.1.2](https://github.com/CAST-genomics/haptools/compare/v0.1.1...v0.1.2) (2023-02-02)


### Bug Fixes

* add poetry readme to fix long_description when publishing to pypi ([#177](https://github.com/CAST-genomics/haptools/issues/177)) ([4050251](https://github.com/CAST-genomics/haptools/commit/405025161ef99b7c18c6adddb58d74bbcff93570))
* checkout in comment bot workflow for CD pipeline ([#181](https://github.com/CAST-genomics/haptools/issues/181)) ([9782d3d](https://github.com/CAST-genomics/haptools/commit/9782d3d5826a7029a74fb845dbfff29a7a0e8f2d))
* checkout ref in comment bot from CD pipeline ([#182](https://github.com/CAST-genomics/haptools/issues/182)) ([e3b92f6](https://github.com/CAST-genomics/haptools/commit/e3b92f6ba2764819d9936c7e37418bcd68469dfb))
* comment bot in CD pipeline ([#180](https://github.com/CAST-genomics/haptools/issues/180)) ([a2f66bd](https://github.com/CAST-genomics/haptools/commit/a2f66bd74dfb1fae8a81614912cba35366fe003b))
* indentation in checkout build ([#183](https://github.com/CAST-genomics/haptools/issues/183)) ([de431ab](https://github.com/CAST-genomics/haptools/commit/de431ab65a15f4dda2e7b6319608331b73ee5771))

## [0.1.1](https://github.com/CAST-genomics/haptools/compare/v0.1.0...v0.1.1) (2023-02-02)


### Documentation

* add dependency to conda env build ([#175](https://github.com/CAST-genomics/haptools/issues/175)) ([6626049](https://github.com/CAST-genomics/haptools/commit/66260498cbaf7fb652ababd6241cdfbc8a865dfd))

## [0.1.0](https://github.com/CAST-genomics/haptools/compare/v0.0.3...v0.1.0) (2023-01-31)


### Features

* a new `--no-normalize` parameter for `simphenotype` ([#156](https://github.com/CAST-genomics/haptools/issues/156)) ([24a0867](https://github.com/CAST-genomics/haptools/commit/24a0867ee8c4b1f495f5dde43600758daa2d916d))
* add `--seed` to `simphenotype` ([#162](https://github.com/CAST-genomics/haptools/issues/162)) ([9f7890a](https://github.com/CAST-genomics/haptools/commit/9f7890a6efd5be90717aaeeeab1f638512c6e2c4))
* change default verbosity to INFO ([#157](https://github.com/CAST-genomics/haptools/issues/157)) ([18ff090](https://github.com/CAST-genomics/haptools/commit/18ff090b431891b16b4aae9d297983d481a3d248))
* github actions for publishing ([#141](https://github.com/CAST-genomics/haptools/issues/141)) ([5cf6f29](https://github.com/CAST-genomics/haptools/commit/5cf6f29b64a49f9405ab78eebab5e0de20963376))
* Significantly decreased runtime of simgenotype  ([#163](https://github.com/CAST-genomics/haptools/issues/163)) ([de011d0](https://github.com/CAST-genomics/haptools/commit/de011d0106a4b50e0ea1e38eff6c567174a3a72f))
* updated `.vcf.gz` output in simgenotype ([#150](https://github.com/CAST-genomics/haptools/issues/150)) ([f61f613](https://github.com/CAST-genomics/haptools/commit/f61f613c55638efb610fb9ca03b0292f42dc2623))
* Updated simgenotype to allow only breakpoint generation ([#167](https://github.com/CAST-genomics/haptools/issues/167)) ([c0c3c97](https://github.com/CAST-genomics/haptools/commit/c0c3c97384c94056f79786f5297251bf9c29f2b9))


### Bug Fixes

* Added an additional test file containing only autosome centromeres ([#168](https://github.com/CAST-genomics/haptools/issues/168)) ([ca220e7](https://github.com/CAST-genomics/haptools/commit/ca220e7319f4b4d6594bad0042a9eb35c2a2d607))
* Added utilization of logging class to karyogram and simgenotype ([#159](https://github.com/CAST-genomics/haptools/issues/159)) ([bcf778a](https://github.com/CAST-genomics/haptools/commit/bcf778a1180d6208c7f24a6f63af1df830734dc5))
* case where genotypes are all the same but heritability is specified in `simphenotype` ([#145](https://github.com/CAST-genomics/haptools/issues/145)) ([063b411](https://github.com/CAST-genomics/haptools/commit/063b411302e4f299f02de82ad4040d2eb3fa8a61))
* Fixed logic error in simgenotype and docs errors ([#172](https://github.com/CAST-genomics/haptools/issues/172)) ([18dd27f](https://github.com/CAST-genomics/haptools/commit/18dd27fbe45c4af0081aa36928796cd61e30041e))
* issue where breakpoints weren't outputting ([#144](https://github.com/CAST-genomics/haptools/issues/144)) ([ff90bff](https://github.com/CAST-genomics/haptools/commit/ff90bff662b5b3ebc98278cf64d4919fbe4699c5))
* logging so that it doesn't affect the root logger ([#154](https://github.com/CAST-genomics/haptools/issues/154)) ([c1dc1ef](https://github.com/CAST-genomics/haptools/commit/c1dc1efed132f365c9268978be37bf25b4ef10b5))
* pin numpy to ensure a recent version of numpy ([#151](https://github.com/CAST-genomics/haptools/issues/151)) ([7c84a45](https://github.com/CAST-genomics/haptools/commit/7c84a454dfe3b57e559a3f4439b1f2dc390700d0))
* Update maps.rst ([#170](https://github.com/CAST-genomics/haptools/issues/170)) ([8ace893](https://github.com/CAST-genomics/haptools/commit/8ace8933644e98eda9483d47f4a895fe8cdd7831))


### Documentation

* better showcase example files and their locations in `simgenotype` ([#166](https://github.com/CAST-genomics/haptools/issues/166)) ([6995407](https://github.com/CAST-genomics/haptools/commit/69954075e547917882f72ca70776169254d16b99))
* updates for v0.1.0 ([#171](https://github.com/CAST-genomics/haptools/issues/171)) ([714142f](https://github.com/CAST-genomics/haptools/commit/714142f06204b284571f5b01c8a2bed3892dada7))
