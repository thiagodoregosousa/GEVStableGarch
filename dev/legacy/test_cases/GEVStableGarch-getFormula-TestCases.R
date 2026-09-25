
# Copyrights (C) 2014 Thiago do Rego Sousa <thiagoestatistico@gmail.com>

# This library is free software; you can redistribute it and/or
# modify it under the terms of the GNU Library General Public
# License as published by the Free Software Foundation; either
# version 2 of the License, or (at your option) any later version.
#
# This library is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Library General Public License for more details.
#
# You should have received a copy of the GNU Library General
# Public License along with this library; if not, write to the
# Free Foundation, Inc., 59 Temple Place, Suite 330, Boston,
# MA  02111-1307  USA



################################################################################
# TEST CASES FOF FUNCTION:               SPECIFICATION:
#  .getFormula                Tests on different types of formulas as input. 
#                                   				                							               
################################################################################


# Test Cases for function 

# expects: formula.mean = ~arma(1,1); formula.variance = ~garch(2,2)
.getFormula(~arma(1,1)+garch(2,2)) 

# expects: formula.mean = ~arma(5,5); formula.variance = ~garch(1,1)
.getFormula(~arma(5,5)+aparch(1,1))

# expects: formula.mean = ~arma(0,0); formula.variance = ~aparch(1,1)
.getFormula(~aparch(1,1))

# expects: formula.mean = ~arma(0,0); formula.variance = ~garch(1,1)
.getFormula(~garch(1,1))

# expects: formula.mean = ~arma(1,1); formula.variance = ~garch(0,0)
.getFormula(~arma(1,1))

# expects: formula.mean = ~arma(1,1); formula.variance = ~garch(0,0)
.getFormula(~arma(0,1))

# expects: formula.mean = ~arma(1,1); formula.variance = ~garch(0,0)
.getFormula(~arma(1,0))

# expects: error
.getFormula(~ar(1))

# expects: error
.getFormula(~ma(1))

# expects: error
.getFormula(~arch(1))

# expects: error
.getFormula(~garch(1,1)+arma(1,1))

# expects: error
.getFormula(~arma(5,5)+aparch(1,1)+arma(2,2))

# expects: error
.getFormula(~aparch(0,1))

# expects: error
.getFormula(~garch(0,1))
####################################################################


