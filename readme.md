# C++ Ported Version of .Net Core 3.1 Financial

C++ financial formulas in single header from most trustworthy source code [Financial.vb](https://github.com/dotnet/corefx/blob/release/3.1/src/Microsoft.VisualBasic.Core/src/Microsoft/VisualBasic/Financial.vb) as part of .Net Core 3.1. 


## Formulas Inside

* DDB

* FV

* IPmt

* IRR

* MIRR

* NPer

* NPV

* Pmt

* PPmt

* PV

* Rate

* SLN

* SYD

## usage

```

#include "include/financial.h"

using namespace  MS_VB_Financial;

double irr = Financial::IRR(yourVectorOfCashflow, guess);

```

## To Do

1. cross platform build

2. add more test cased from [FinancialTests.cs](https://github.com/dotnet/corefx/blob/release/3.1/src/Microsoft.VisualBasic.Core/tests/FinancialTests.cs)

3. add performance test

4. build nodejs addon

5. add github action