# CASM-Smart-Phase

| Master                                        | Develop                                         |
| --------------------------------------------- | ----------------------------------------------- |
| [![Main Badge][circle-master]][circle-master-base] | [![Develop Badge][circle-develop]][circle-develop-base] |

- [Smart-Phase](#smart-phase)
  - [Running Smart-Phase](#running-smart-phase)
- [Python utility scripts](#python-utility-scripts)
- [Docker image](#docker-image)
- [LICENCE](#licence)

## Smart-Phase

The docker image also contains [Smart-Phase]
Please note, we do not support the Smart-Phase software, please contact the Smart-Phase repository owner.

### Running Smart-Phase

In order to run smart-phase in the docker image either use the default command:

```bash
docker run casm-smart-phase:latest

usage: Welcome to SmartPhase! A dedicated tool designed to assist in the
       rapid and accurate phasing of variant combinations for clinical
       analysis. Please refer to the following list of options on how to
       pass the necessary parameters for use:
       ....
```

Or

```bash
docker run smart-phase:latest java -jar /opt/wsi-t78/smartPhase.jar

usage: Welcome to SmartPhase! A dedicated tool designed to assist in the
       rapid and accurate phasing of variant combinations for clinical
       analysis. Please refer to the following list of options on how to
       pass the necessary parameters for use:
       ....
```

## Python utility scripts

The [python](/python) subdirectory contains utility scripts to merge [Smart-Phase] MNVs
For more detail see the section [README](/python/README.md).

## Docker image

Eventually this base directory will contain the Docker script for a [Smart-Phase] and
utility script container.

## LICENCE

```txt
Copyright (c) 2021 Genome Research Ltd.

Author: CASM/Cancer IT <cgphelp@sanger.ac.uk>

This file is part of CASM-Smart-Phase.

This program is free software: you can redistribute it and/or modify it under
the terms of the GNU Affero General Public License as published by the Free
Software Foundation; either version 3 of the License, or (at your option) any
later version.

This program is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
FOR A PARTICULAR PURPOSE. See the GNU Affero General Public License for more
details.

You should have received a copy of the GNU Affero General Public License
along with this program. If not, see <http://www.gnu.org/licenses/>.

1. The usage of a range of years within a copyright statement contained within
this distribution should be interpreted as being equivalent to a list of years
including the first and last year specified and all consecutive years between
them. For example, a copyright statement that reads ‘Copyright (c) 2005, 2007-
2009, 2011-2012’ should be interpreted as being identical to a statement that
reads ‘Copyright (c) 2005, 2007, 2008, 2009, 2011, 2012’ and a copyright
statement that reads ‘Copyright (c) 2005-2012’ should be interpreted as being
identical to a statement that reads ‘Copyright (c) 2005, 2006, 2007, 2008,
2009, 2010, 2011, 2012’."
```

<!-- Reference Links -->

[circle-develop]: https://circleci.com/gh/cancerit/CASM-Smart-Phase/tree/develop.svg?style=shield
[circle-develop-base]: https://circleci.com/gh/cancerit/CASM-Smart-Phase/tree/develop
[circle-master]: https://circleci.com/gh/cancerit/CASM-Smart-Phase/tree/main.svg?style=shield
[circle-master-base]: https://circleci.com/gh/cancerit/CASM-Smart-Phase/tree/main
[smart-phase]: https://github.com/paulhager/smart-phase
