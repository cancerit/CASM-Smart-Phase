#!/usr/bin/env bash

##   LICENSE
# Copyright (c) 2018-2019 Genome Research Ltd.
# Author: Cancer Genome Project cgphelp@sanger.ac.uk
#
#
# This file is part of vcf_flag_modifier.
#
# vcf_flag_modifier is free software: you can redistribute it and/or modify it under
# the terms of the GNU Affero General Public License as published by the Free
# Software Foundation; either version 3 of the License, or (at your option) any
# later version.
#
# This program is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
# FOR A PARTICULAR PURPOSE. See the GNU Affero General Public License for more
# details.
#
# You should have received a copy of the GNU Affero General Public License
# along with this program. If not, see <http://www.gnu.org/licenses/>.
#
#    1. The usage of a range of years within a copyright statement contained within
#    this distribution should be interpreted as being equivalent to a list of years
#    including the first and last year specified and all consecutive years between
#    them. For example, a copyright statement that reads ‘Copyright (c) 2005, 2007-
#    2009, 2011-2012’ should be interpreted as being identical to a statement that
#    reads ‘Copyright (c) 2005, 2007, 2008, 2009, 2011, 2012’ and a copyright
#    statement that reads ‘Copyright (c) 2005-2012’ should be interpreted as being
#    identical to a statement that reads ‘Copyright (c) 2005, 2006, 2007, 2008,
#    2009, 2010, 2011, 2012’."
#
#

set -e

SRC_DIR=`dirname $0`/casmsmartphase
TEST_DIR=`dirname $0`/tests

 #!/usr/bin/env bash
set -e
pytest $TEST_DIR \
 --cov-branch\
 --cov-report term\
 --cov-report html\
 --cov-fail-under=89\
 --cov=casmsmartphase\
 --cov=scripts\
 -x
set +e

# these should not die:

echo -e "\n#################################"
echo      "# Running pycodestyle (style)   #"
echo      "#################################"
pycodestyle casmsmartphase

echo -e "\n#########################################"
echo      "# Running radon (cyclomatic complexity) #"
echo      "#########################################"
radon cc -nc casmsmartphase

echo -e "\n#########################################"
echo      "# Running radon (maintainability index) #"
echo      "#########################################"
radon mi -s casmsmartphase

echo -e "\n##############################"
echo      "# Running mdl (markdownlint) #"
echo      "##############################"
mdl .

exit 0 # don't die based on assements of code quality
