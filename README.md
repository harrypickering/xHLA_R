xHLA: Fast and accurate HLA typing from short read sequence data
================================================================

Installation
------------
Compile docker image:

```bash
$ cd docker
$ make build
$ make deploy
```

Usage
-----
Run xHLA caller directly on a BAM file.

```bash
docker run -v `pwd`:`pwd` -w `pwd` docker-dev.hli.io/xchao/hla \
    --sample_id test --input_bam_path test.bam \
    --output_path test
```

License
-------
The HLA Typing Software Code (the "Code") is made available by Human
Longevity, Inc. ("HLI") on a non-exclusive, non-sublicensable,
non-transferable basis solely for non-commercial academic research use.
Commercial use of the Code is expressly prohibited.  If you would like to obtain
a license to the Code for commercial use, please contact HLI at
bizdev@humanlongevity.com.  HLI MAKES NO REPRESENTATIONS OR WARRANTIES
WHATSOEVER, EITHER EXPRESS OR IMPLIED, WITH RESPECT TO THE CODE PROVIDED
HEREUNDER. IMPLIED WARRANTIES OF MERCHANTABILITY OR FITNESS FOR A PARTICULAR
PURPOSE WITH RESPECT TO CODE ARE EXPRESSLY DISCLAIMED. THE CODE IS FURNISHED
"AS IS" AND "WITH ALL FAULTS" AND DOWNLOADING OR USING THE CODE
IS UNDERTAKEN AT YOUR OWN RISK.  TO THE FULLEST EXTENT ALLOWED BY APPLICABLE
LAW, IN NO EVENT SHALL HLI BE LIABLE, WHETHER IN CONTRACT, TORT, WARRANTY, OR
UNDER ANY STATUTE OR ON ANY OTHER BASIS FOR SPECIAL, INCIDENTAL, INDIRECT,
PUNITIVE, MULTIPLE OR CONSEQUENTIAL DAMAGES SUSTAINED BY YOU OR ANY OTHER PERSON
OR ENTITY ON ACCOUNT OF USE OR POSSESSION OF THE CODE, WHETHER OR NOT
FORESEEABLE AND WHETHER OR NOT HLI HAS BEEN ADVISED OF THE POSSIBILITY OF SUCH
DAMAGES, INCLUDING WITHOUT LIMITATION DAMAGES ARISING FROM OR RELATED TO LOSS OF
USE, LOSS OF DATA, DOWNTIME, OR FOR LOSS OF REVENUE, PROFITS, GOODWILL, BUSINESS
OR OTHER FINANCIAL LOSS.