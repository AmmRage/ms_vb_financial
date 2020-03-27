#pragma once

#include <vector>
#include <cmath>
#include <stdexcept>
#include "exceptionhelper.h"

// Copyright (c) Microsoft Corporation.  All rights reserved.



namespace MS_VB_Financial
{
	///
	///Indicates when payments are due when calling financial methods.
	///指示在调用财务方法时应何时付款。
	///
	enum class DueDate {
		///
		///Falls at the beginning of the date interval
		///
	    BegOfPeriod = 1,

	    ///
	    ///Falls at the end of the date interval
	    ///
		EndOfPeriod = 0,
	};

	class Financial
	{

		// ============================================================================
		// Financial functions.
		// ============================================================================

	private:
		static constexpr double cnL_IT_STEP = 0.00001;
		static constexpr double cnL_IT_EPSILON = 0.0000001;


		// -------------------------------------------------------------
		// 
		// Name          : DDB

		// Purpose       : Calculates depreciation for a period based on the
		// doubly-declining balance method. Returns the result.
		// It raises an error if parameters are invalid.

		// Derivation    : At each period, 2*balance/nper is subtracted
		// from the balance.  The balance starts at the purchase
		// value, and the total of the payments may not
		// exceed (value - salvage).  The algorithm uses a non-
		// iterative method to calculate the payment.
		// Note that only the integral values of nper and per make any
		// sense.  However, Excel allowed non-integral values, and
		// thus these routines also work with non-integral input.
		// 
		// PMT = 2 * (value / nper) * ( (nper -2) / nper ) ^ (per - 1)
		// 
		// total = value * ( 1 - ( (nper - 2) / nper ) ^ per )
		// 
		// excess = total - (value - salvage)
		// 
		// ddb =  PMT        : if  excess <= 0
		// PMT-excess : if  PMT >= excess > 0
		// 0          : if  excess > PMT
		// Returns       : Double
		// 
		// -------------------------------------------------------------
		// 
	public:
		static double DDB(double Cost, double Salvage, double Life, double Period, double Factor = 2.0)
		{
			double dRet;
			double dTot;
			double dExcess;
			double dTemp;
			double dNTemp;

			// Handle invalid parameters
			if (Factor <= 0.0 || Salvage < 0.0 || Period <= 0.0 || Period > Life)
			{
				throw std::invalid_argument(("Argument_InvalidValue1: Factor"));
			}

			// Handle special (trivial) cases
			if (Cost <= 0.0)
			{
				return 0.0;
			}

			if (Life < 2.0)
			{
				return (Cost - Salvage);
			}

			if (Life == 2.0 && Period > 1.0)
			{
				return 0.0;
			}

			if (Life == 2.0 && Period <= 1.0)
			{
				return (Cost - Salvage);
			}

			if (Period <= 1.0)
			{
				dRet = Cost * Factor / Life;
				dTemp = Cost - Salvage;
				if (dRet > dTemp)
				{
					return dTemp;
				}
				else
				{
					return dRet;
				}
			}

			// Perform the calculation
			dTemp = (Life - Factor) / Life;
			dNTemp = Period - 1.0;

			// WARSI Using the exponent operator for pow(..) in C code of DDB. Still got
			// to make sure that they (pow and ^) are same for all conditions
			dRet = Factor * Cost / Life * std::pow(dTemp, dNTemp);

			// WARSI Using the exponent operator for pow(..) in C code of DDB. Still got
			// to make sure that they (pow and ^) are same for all conditions
			dTot = Cost * (1 - std::pow(dTemp, Period));
			dExcess = dTot - Cost + Salvage;

			if (dExcess > 0.0)
			{
				dRet = dRet - dExcess;
			}

			if (dRet >= 0.0)
			{
				return dRet;
			}
			else
			{
				return 0.0;
			}
		}



		// -------------------------------------------------------------
		// 
		// Name                      : FV
		// Purpose                   : It is computed as -

		// (1+rate)^nper - 1
		// fv = -pv*(1+rate)^nper - PMT*(1+rate*type)* -----------------
		// rate
		// 
		// fv = -pv - PMT * nper        : if rate == 0
		// 
		// 
		// Returns                   : Double
		// 
		// -------------------------------------------------------------
		// 
		static double FV(double Rate, double NPer, double Pmt, double PV = 0, DueDate Due = DueDate::EndOfPeriod)
		{
			return FV_Internal(Rate, NPer, Pmt, PV, Due);
		}



	private:
		static double FV_Internal(double Rate, double NPer, double Pmt, double PV = 0, DueDate Due = DueDate::EndOfPeriod)
		{
			double dTemp;
			double dTemp2;
			double dTemp3;

			// Performing calculation
			if (Rate == 0.0)
			{
				return (-PV - Pmt * NPer);
			}

			if (Due != DueDate::EndOfPeriod)
			{
				dTemp = 1.0 + Rate;
			}
			else
			{
				dTemp = 1.0;
			}

			dTemp3 = 1.0 + Rate;
			dTemp2 = std::pow(dTemp3, NPer);

			// Do divides before multiplies to avoid OverflowExceptions
			return ((-PV) * dTemp2) - ((Pmt / Rate) * dTemp * (dTemp2 - 1.0));
		}



		// -------------------------------------------------------------
		// 
		// Name       : IPmt
		// Purpose    : This function calculates the interest part of a
		// payment for a given period.  The payment is part of
		// a series of regular payments described by the other
		// arguments.  The value is returned.  The function
		// Raises an expection if params are invalid. This function
		// calls FV and PMT. It calculates value of annuity (FV) at the
		// begining of period for which IPMT is desired. Since interest
		// rate is constant FV*rate would give IPMT.
		// 
		// if type = 1 and per = 1 then IPMT = 0.
		// 
		// if (type = 0 ) IPMT = FV(per-1)*rate
		// else IPMT = FV(per-2)*rate
		// Returns    : Double
		// 
		// -------------------------------------------------------------
		// 
	public:
		static double IPmt(double Rate, double Per, double NPer, double PV, double FV = 0, DueDate Due = DueDate::EndOfPeriod)
		{
			double Pmt;
			double dTFv;
			double dTemp;

			if (Due != DueDate::EndOfPeriod)
			{
				dTemp = 2.0;
			}
			else
			{
				dTemp = 1.0;
			}

			// Type = 0 or non-zero only. Offset to calculate FV
			if ((Per <= 0.0) || (Per >= NPer + 1))
			{
				throw std::invalid_argument(("Argument_InvalidValue1: Per"));
			}

			if ((Due != DueDate::EndOfPeriod) && (Per == 1.0))
			{
				return 0.0;
			}

			// Calculate PMT (i.e. annuity) for given parms. Rqrd for FV
			Pmt = PMT_Internal(Rate, NPer, PV, FV, Due);

			// Calculate FV just before the period on which interest would be applied
			if (Due != DueDate::EndOfPeriod)
			{
				PV = PV + Pmt;
			}

			dTFv = FV_Internal(Rate, (Per - dTemp), Pmt, PV, DueDate::EndOfPeriod);

			return (dTFv * Rate);
		}



		// -------------------------------------------------------------
		// 
		// Name                      : IRR
		// Purpose                   : This function uses an iterative procedure to find
		// the Internal Rate of Return of an investment.  The algorithm
		// basically uses the secant method to find a rate for which the
		// NPV of the cash flow is 0.
		// This function raises an exception if the parameters are invalid.
		// 
		// This routine uses a slightly different version of the secant
		// routine in Rate.  The basic changes are:
		// - uses LDoNPV to get the 'Y-value'
		// - does not allow Rate to go below -1.
		// (if the Rate drops below -1, it is forced above again)
		// - has a double condition for termination:
		// NPV = 0 (within L_IT_EPSILON)
		// Rate1 - Rate0  approaches zero (rate is converging)
		// This last does not parallel Excel, but avoids a class of
		// spurious answers.  Otherwise, performance is comparable to
		// Excel's, and accuracy is often better.
		// 
		// Returns                   : Double
		// 
		// -------------------------------------------------------------
		// 
		static double IRR(std::vector<double> &ValueArray, double Guess = 0.1)
		{
			double dTemp;
			double dRate0;
			double dRate1;
			double dNPv0;
			double dNpv1;
			double dNpvEpsilon;
			double dTemp1;
			int lIndex;
			int lCVal;
			int lUpper;

			// Compiler assures that rank of ValueArray is always 1, no need to check it.  
			// WARSI Check for error codes returned by UBound. Check if they match with C code
			try // Needed to catch dynamic arrays which have not been constructed yet.
			{
				lUpper = ValueArray.size()-1;
			}
			catch (const StackOverflowException &ex)
			{
				throw ex;
			}
			catch (const OutOfMemoryException &ex)
			{
				throw ex;
			}
			catch (const ThreadAbortException &ex)
			{
				throw ex;
			}
			catch (...)
			{
				throw std::invalid_argument("Argument_InvalidValue1");
			}

			lCVal = lUpper + 1;

			// Function fails for invalid parameters
			if (Guess <= (-1.0))
			{
				throw std::invalid_argument(("Argument_InvalidValue1: Guess"));
			}

			if (lCVal <= 1)
			{
				throw std::invalid_argument(("Argument_InvalidValue1: ValueArray"));
			}

			// We scale the epsilon depending on cash flow values. It is necessary
			// because even in max accuracy case where diff is in 16th digit it
			// would get scaled up.
			if (ValueArray[0] > 0.0)
			{
				dTemp = ValueArray[0];
			}
			else
			{
				dTemp = -ValueArray[0];
			}

			for (lIndex = 0; lIndex <= lUpper; lIndex++)
			{
				// Get max of values in cash flow
				if (ValueArray[lIndex] > dTemp)
				{
					dTemp = ValueArray[lIndex];
				}
				else if ((-ValueArray[lIndex]) > dTemp)
				{
					dTemp = -ValueArray[lIndex];
				}
			}

			dNpvEpsilon = dTemp * cnL_IT_EPSILON * 0.01;

			// Set up the initial values for the secant method
			dRate0 = Guess;
			dNPv0 = OptPV2(ValueArray, dRate0);

			if (dNPv0 > 0.0)
			{
				dRate1 = dRate0 + cnL_IT_STEP;
			}
			else
			{
				dRate1 = dRate0 - cnL_IT_STEP;
			}

			if (dRate1 <= (-1.0))
			{
				throw std::invalid_argument(("Argument_InvalidValue1: Rate"));
			}

			dNpv1 = OptPV2(ValueArray, dRate1);

			for (lIndex = 0; lIndex <= 39; lIndex++)
			{
				if (dNpv1 == dNPv0)
				{
					if (dRate1 > dRate0)
					{
						dRate0 = dRate0 - cnL_IT_STEP;
					}
					else
					{
						dRate0 = dRate0 + cnL_IT_STEP;
					}
					dNPv0 = OptPV2(ValueArray, dRate0);
					if (dNpv1 == dNPv0)
					{
					    throw std::invalid_argument(("Argument_InvalidValue1"));
					}
				}

				dRate0 = dRate1 - (dRate1 - dRate0) * dNpv1 / (dNpv1 - dNPv0);

				// Secant method of generating next approximation
				if (dRate0 <= (-1.0))
				{
					dRate0 = (dRate1 - 1.0) * 0.5;
				}

				// Basically give the algorithm a second chance. Helps the
				// algorithm when it starts to diverge to -ve side
				dNPv0 = OptPV2(ValueArray, dRate0);
				if (dRate0 > dRate1)
				{
					dTemp = dRate0 - dRate1;
				}
				else
				{
					dTemp = dRate1 - dRate0;
				}

				if (dNPv0 > 0.0)
				{
					dTemp1 = dNPv0;
				}
				else
				{
					dTemp1 = -dNPv0;
				}

				// Test : npv - > 0 and rate converges
				if (dTemp1 < dNpvEpsilon && dTemp < cnL_IT_EPSILON)
				{
					return dRate0;
				}

				// Exchange the values - store the new values in the 1's
				dTemp = dNPv0;
				dNPv0 = dNpv1;
				dNpv1 = dTemp;
				dTemp = dRate0;
				dRate0 = dRate1;
				dRate1 = dTemp;
			}

			throw std::invalid_argument(("Argument_InvalidValue1: Rate"));
		}



		// -------------------------------------------------------------
		// 
		// Name                      : MIRR
		// Returns                   : Double
		// 
		// -------------------------------------------------------------
		// 
		static double MIRR(std::vector<double> &ValueArray, double FinanceRate, double ReinvestRate)
		{
			double dNpvPos;
			double dNpvNeg;
			double dTemp;
			double dTemp1;
			double dNTemp2;
			int lCVal;
			int lLower;
			int lUpper;


			lLower = 0;
			lUpper = ValueArray.size()-1;
			lCVal = lUpper - lLower + 1;

			if (FinanceRate == -1.0)
			{
				throw std::invalid_argument(("Argument_InvalidValue1: FinanceRate"));
			}

			if (ReinvestRate == -1.0)
			{
				throw std::invalid_argument(("Argument_InvalidValue1: ReinvestRate"));
			}

			if (lCVal <= 1)
			{
				throw std::invalid_argument(("Argument_InvalidValue1: ValueArray"));
			}

			dNpvNeg = LDoNPV(FinanceRate, ValueArray, -1);
			if (dNpvNeg == 0.0)
			{
				throw std::invalid_argument(("DivideByZeroException: Financial_CalcDivByZero"));
			}

			dNpvPos = LDoNPV(ReinvestRate, ValueArray, 1); // npv of +ve values
			dTemp1 = ReinvestRate + 1.0;
			dNTemp2 = lCVal;

			dTemp = -dNpvPos * std::pow(dTemp1, dNTemp2) / (dNpvNeg * (FinanceRate + 1.0));

			if (dTemp < 0.0)
			{
				throw std::invalid_argument(("Argument_InvalidValue1"));
			}

			dTemp1 = 1 / (lCVal - 1.0);

			return std::pow(dTemp, dTemp1) - 1.0;
		}



		// -------------------------------------------------------------
		// 
		// Name                      : NPer
		// Purpose                   :
		// 
		// -fv + PMT*(1+rate*type) / rate
		// (1+rate)^nper = ------------------------------
		// pv + PMT*(1+rate*type) / rate
		// 
		// this yields the log expression used in this function.
		// 
		// nper = (-fv - pv) / PMT     : if rate == 0
		// 
		// Returns                   : Double
		// 
		// -------------------------------------------------------------
		// 
		static double NPer(double Rate, double Pmt, double PV, double FV = 0, DueDate Due = DueDate::EndOfPeriod)
		{
			double dTemp3;
			double dTempFv;
			double dTempPv;
			double dTemp4;

			// Checking Error Conditions
			if (Rate <= -1.0)
			{
				throw std::invalid_argument(("Argument_InvalidValue1: Rate"));
			}

			if (Rate == 0.0)
			{
				if (Pmt == 0.0)
				{
				throw std::invalid_argument(("Argument_InvalidValue1: Pmt"));
				}
				return (-(PV + FV) / Pmt);
			}
			else
			{
				if (Due != DueDate::EndOfPeriod)
				{
					dTemp3 = Pmt * (1.0 + Rate) / Rate;
				}
				else
				{
					dTemp3 = Pmt / Rate;
				}
				dTempFv = -FV + dTemp3;
				dTempPv = PV + dTemp3;

				// Make sure the values fit the domain of log()
				if (dTempFv < 0.0 && dTempPv < 0.0)
				{
					dTempFv = -1 * dTempFv;
					dTempPv = -1 * dTempPv;
				}
				else if (dTempFv <= 0.0 || dTempPv <= 0.0)
				{
				throw std::invalid_argument(("Argument_InvalidValue1: Financial_CannotCalculateNPer"));
				}

				dTemp4 = Rate + 1.0;
				return (std::log(dTempFv) - std::log(dTempPv)) / std::log(dTemp4);
			}
		}



		// -------------------------------------------------------------
		// 
		// Name       : NPV
		// Purpose                   :
		// This function calculates the Net Present Value of a series of
		// payments at a given rate.  It uses LDoNPV to get the value.  No
		// real work is done here, just some error checking.
		// As with the others, this function puts its status in *lpwStatus,
		// and returns the result as a double.
		// 
		// Value1       Value2      Value3
		// npv = -------- + ---------- + ---------- + ... for the series...
		// (1+rate)   (1+rate)^2   (1+rate)^3
		// 
		// 
		// Returns    : Double
		// 
		// -------------------------------------------------------------
		// 
		static double NPV(double Rate, std::vector<double> &ValueArray)
		{
			int lCVal;
			int lLower;
			int lUpper;

			if ((ValueArray.empty()))
			{
				throw std::invalid_argument(("Argument_InvalidValue1: ValueArray"));
			}

			lLower = 0;
			lUpper = ValueArray.size()-1;
			lCVal = lUpper - lLower + 1;

			if (Rate == (-1.0))
			{
				throw std::invalid_argument(("Argument_InvalidValue1: Rate"));
			}
			if (lCVal < 1)
			{
				throw std::invalid_argument(("Argument_InvalidValue1: ValueArray"));
			}

			return LDoNPV(Rate, ValueArray, 0);
		}



		// -------------------------------------------------------------
		// 
		// Name                      : PMT
		// Purpose                   :
		// This function, together with the four following
		// it (Pv, Fv, NPer and Rate), can calculate
		// a certain value associated with a regular series of
		// equal.size()-1d payments.  This series can be fully described
		// by these values:
		// Pv   - present value
		// Fv   - future value (at end of series)
		// PMT  - the regular payment
		// nPer - the number of 'periods' over which the
		// money is paid
		// Rate - the interest rate per period.
		// (type - payments at beginning (1) or end (0) of
		// the period).
		// Each function can determine one of the values, given the others.
		// 
		// General Function for the above values:
		// 
		// (1+rate)^nper - 1
		// pv * (1+rate)^nper + PMT*(1+rate*type)*----------------- + fv  = 0
		// rate
		// rate == 0  ->  pv + PMT*nper + fv = 0
		// 
		// Thus:
		// (-fv - pv*(1+rate)^nper) * rate
		// PMT = -------------------------------------
		// (1+rate*type) * ( (1+rate)^nper - 1 )
		// 
		// PMT = (-fv - pv) / nper    : if rate == 0
		// 
		// 
		// Returns                   : Double
		// 
		// -------------------------------------------------------------
		// 
		static double Pmt(double Rate, double NPer, double PV, double FV = 0, DueDate Due = DueDate::EndOfPeriod)
		{
			return PMT_Internal(Rate, NPer, PV, FV, Due);
		}

	private:
		static double PMT_Internal(double Rate, double NPer, double PV, double FV = 0, DueDate Due = DueDate::EndOfPeriod)
		{
			double dTemp;
			double dTemp2;
			double dTemp3;

			// Checking for error conditions
			if (NPer == 0.0)
			{
				throw std::invalid_argument(("Argument_InvalidValue1: NPer"));
			}

			if (Rate == 0.0)
			{
				return ((-FV - PV) / NPer);
			}
			else
			{
				if (Due != DueDate::EndOfPeriod)
				{
					dTemp = 1.0 + Rate;
				}
				else
				{
					dTemp = 1.0;
				}
				dTemp3 = Rate + 1.0;
				// WARSI Using the exponent operator for pow(..) in C code of PMT. Still got
				// to make sure that they (pow and ^) are same for all conditions
				dTemp2 = std::pow(dTemp3, NPer);
				return ((-FV - PV * dTemp2) / (dTemp * (dTemp2 - 1.0)) * Rate);
			}
		}



		// -------------------------------------------------------------
		// 
		// Name                      : PPmt
		// Purpose                   : This function calculates the principal part of a
		// payment for a given period.
		// 
		// Since PMT = IPMT +PPMT therefore
		// PPMT = PMT - IPMT
		// 
		// Returns                   : Double
		// 
		// -------------------------------------------------------------
		// 
	public:
		static double PPmt(double Rate, double Per, double NPer, double PV, double FV = 0, DueDate Due = DueDate::EndOfPeriod)
		{
			double Pmt;
			double dIPMT;

			// Checking for error conditions
			if ((Per <= 0.0) || (Per >= (NPer + 1)))
			{
				throw std::invalid_argument(("Argument_InvalidValue1: Per"));
			}

			Pmt = PMT_Internal(Rate, NPer, PV, FV, Due);
			dIPMT = IPmt(Rate, Per, NPer, PV, FV, Due);

			return (Pmt - dIPMT);
		}



		// -------------------------------------------------------------
		// 
		// Name       : PV
		// Purpose    :
		// -fv - PMT * (1+rate*type) * ( (1+rate)^nper-1) / rate
		// pv = -----------------------------------------------------
		// (1 + rate) ^ nper
		// 
		// pv = -fv - PMT * nper     : if rate == 0
		// Returns    : Double
		// 
		// -------------------------------------------------------------
		// 
		static double PV(double Rate, double NPer, double Pmt, double FV = 0, DueDate Due = DueDate::EndOfPeriod)
		{
			double dTemp;
			double dTemp2;
			double dTemp3;

			if (Rate == 0.0)
			{
				return (-FV - Pmt * NPer);
			}
			else
			{
				if (Due != DueDate::EndOfPeriod)
				{
					dTemp = 1.0 + Rate;
				}
				else
				{
					dTemp = 1.0;
				}
				dTemp3 = 1.0 + Rate;

				// WARSI Using the exponent operator for pow(..) in C code of PV. Still got
				// to make sure that they (pow and ^) are same for all conditions
				dTemp2 = std::pow(dTemp3, NPer);

				// Do divides before multiplies to avoid OverFlowExceptions
				return (-(FV + Pmt * dTemp * ((dTemp2 - 1.0) / Rate)) / dTemp2);
			}
		}



		// -------------------------------------------------------------
		// 
		// Name       : Rate
		// Purpose    :
		// See PMT, above, for general details.  This
		// function is not as simple as the others.  Due to the
		// nature of the equation that links the 5 values (see
		// Excel manual - PV), it is not practical to solve for
		// rate algebraically.  As a result, this function implements
		// the secant method of approximation.  LEvalRate provides
		// the 'Y-values', for given rates.
		// 
		// Basic secant method:
		// determine Rate0 and Rate1.  Use LEvalRate to get Y0, Y1.
		// 
		// Y0
		// Rate2 = Rate1 - (Rate1 - Rate0) * ---------
		// (Y1 - Y0)
		// 
		// Get Y2 from Rate2, LEvalRate.  move 1->0, 2->1 and repeat.
		// 
		// stop when abs( Yn ) < L_IT_EPSILON
		// 
		// 
		// Returns    : Double
		// 
		// -------------------------------------------------------------
		// 
		static double Rate(double NPer, double Pmt, double PV, double FV = 0, DueDate Due = DueDate::EndOfPeriod, double Guess = 0.1)
		{
			double dTemp;
			double dRate0;
			double dRate1;
			double dY0;
			double dY1;
			int I;

			// Check for error condition
			if (NPer <= 0.0)
			{
				throw std::invalid_argument(("Argument_InvalidValue1: Rate_NPerMustBeGTZero"));
			}

			dRate0 = Guess;
			dY0 = LEvalRate(dRate0, NPer, Pmt, PV, FV, Due);
			if (dY0 > 0)
			{
				dRate1 = (dRate0 / 2);
			}
			else
			{
				dRate1 = (dRate0 * 2);
			}

			dY1 = LEvalRate(dRate1, NPer, Pmt, PV, FV, Due);

			for (I = 0; I <= 39; I++)
			{
				if (dY1 == dY0)
				{
					if (dRate1 > dRate0)
					{
						dRate0 = dRate0 - cnL_IT_STEP;
					}
					else
					{
						dRate0 = dRate0 - cnL_IT_STEP * (-1);
					}
					dY0 = LEvalRate(dRate0, NPer, Pmt, PV, FV, Due);
					if (dY1 == dY0)
					{
					    throw std::invalid_argument(("Financial_CalcDivByZero"));
					}
				}

				dRate0 = dRate1 - (dRate1 - dRate0) * dY1 / (dY1 - dY0);

				// Secant method of generating next approximation
				dY0 = LEvalRate(dRate0, NPer, Pmt, PV, FV, Due);
				if (std::abs(dY0) < cnL_IT_EPSILON)
				{
					return dRate0;
				}

				dTemp = dY0;
				dY0 = dY1;
				dY1 = dTemp;
				dTemp = dRate0;
				dRate0 = dRate1;
				dRate1 = dTemp;
			}

				throw std::invalid_argument(("Argument_InvalidValue1: Financial_CannotCalculateRate"));
		}



		// -------------------------------------------------------------
		// 
		// Name                      : SLN
		// Purpose                   : It calculates the depreciation by the straight
		// line method and returns the result. It raises
		// an error if parameters are invalid.
		// 
		// sln = (value - salvage) / nper
		// 
		// Returns                   : Double
		// 
		// -------------------------------------------------------------
		// 
		static double SLN(double Cost, double Salvage, double Life)
		{
			if (Life == 0.0)
			{
				throw std::invalid_argument(("Argument_InvalidValue1: Financial_LifeNEZero"));
			}

			return (Cost - Salvage) / (Life);
		}



		// -------------------------------------------------------------
		// 
		// Name       : SYD
		// Purpose    : Calculates depreciation for a period by the
		// sum-of-years-digits method. The result is returned.
		// It raises an error if parameters are invalid.
		// 
		// 2
		// syd = (value - salvage) (nper - per + 1) * ------------
		// (nper)(nper + 1)
		// 
		// Derivation : The value of the asset is divided into even parts.
		// The first period gets N, the second gets N-1, the last
		// gets 1.
		// Returns    : Double
		// 
		// -------------------------------------------------------------
		// 
		static double SYD(double Cost, double Salvage, double Life, double Period)
		{
			double Result;

			if (Salvage < 0.0)
			{
				throw std::invalid_argument(("Argument_InvalidValue1: Salvage"));
			}
			if (Period > Life)
			{
				throw std::invalid_argument(("Argument_InvalidValue1: Financial_PeriodLELife"));
			}
			if (Period <= 0.0)
			{
				throw std::invalid_argument(("Argument_InvalidValue1: Period"));
			}

			// Avoid OverflowExceptions by dividing before multiplying
			Result = (Cost - Salvage) / (Life * (Life + 1));
			return (Result * (Life + 1 - Period) * 2);
		}


		// -------------------------------------------------------------
		// 
		// Name       : LEvalRate
		// Purpose    : A local helper function. Does a useful calculation
		// for Rate.  The function is derived from the General
		// formulation noted above (PMT).
		// Returns    : Double
		// 
		// -------------------------------------------------------------
		// 
	private:
		static double LEvalRate(double Rate, double NPer, double Pmt, double PV, double dFv, DueDate Due)
		{
			double dTemp1;
			double dTemp2;
			double dTemp3;

			if (Rate == 0.0)
			{
				return (PV + Pmt * NPer + dFv);
			}
			else
			{
				dTemp3 = Rate + 1.0;
				// WARSI Using the exponent operator for pow(..) in C code of LEvalRate. Still got
				// to make sure that they (pow and ^) are same for all conditions
				dTemp1 = std::pow(dTemp3, NPer);

				if (Due != DueDate::EndOfPeriod)
				{
					dTemp2 = 1 + Rate;
				}
				else
				{
					dTemp2 = 1.0;
				}
				return (PV * dTemp1 + Pmt * dTemp2 * (dTemp1 - 1) / Rate + dFv);
			}
		}


		// -------------------------------------------------------------
		// 
		// Name                      : LDoNPV
		// Purpose                   :
		// This function performs NPV calculations for NPV,
		// MIRR, and IRR.  The wNType variable is used to set
		// the type of calculation: 0 -> do all values,
		// 1 -> only values > 0,
		// -1 -> only values < 0.
		// Note the array pointer, lpdblVal, is preceded by the count of
		// the entries.
		// Since this is just an internal-use function, no fancy exports
		// are done.  It assumes that error checking is done by the caller.
		// 
		// Value1      Value2       Value3
		// npv = -------- + ---------- + ---------- + ... for the series...
		// (1+rate)   (1+rate)^2   (1+rate)^3
		// 
		// Returns                   : Double
		// 
		// -------------------------------------------------------------
		// 

		static double LDoNPV(double Rate, std::vector<double> &ValueArray, int iWNType)
		{
			bool bSkipPos;
			bool bSkipNeg;

			double dTemp2;
			double dTotal;
			double dTVal;
			int I;
			int lLower;
			int lUpper;

			bSkipPos = iWNType < 0;
			bSkipNeg = iWNType > 0;

			dTemp2 = 1.0;
			dTotal = 0.0;

			lLower = 0;
			lUpper = ValueArray.size()-1;

			for (I = lLower; I <= lUpper; I++)
			{
				dTVal = ValueArray[I];
				dTemp2 = dTemp2 + dTemp2 * Rate;

				if (!((bSkipPos && dTVal > 0.0) || (bSkipNeg && dTVal < 0.0)))
				{
					dTotal = dTotal + dTVal / dTemp2;
				}
			}

			return dTotal;
		}

		// ------------------------------------------------------------------------------------------------------
		// Optimized version of PV2
		// ------------------------------------------------------------------------------------------------------

		static double OptPV2(std::vector<double> &ValueArray, double Guess = 0.1)
		{
			int lUpper, lLower, lIndex;

			lLower = 0;
			lUpper = ValueArray.size()-1;

			double dTotal = 0.0;
			double divRate = 1.0 + Guess;

			while (lLower <= lUpper && ValueArray[lLower] == 0.0)
			{
				lLower = lLower + 1;
			}

			for (lIndex = lUpper; lIndex >= lLower; lIndex += -1)
			{
				dTotal = dTotal / divRate;
				dTotal = dTotal + ValueArray[lIndex];
			}
			return dTotal;
		}
	};
}
