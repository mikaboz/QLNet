﻿//  Copyright (C) 2008-2016 Andrea Maggiulli (a.maggiulli@gmail.com)
//  
//  This file is part of QLNet Project https://github.com/amaggiulli/qlnet
//  QLNet is free software: you can redistribute it and/or modify it
//  under the terms of the QLNet license.  You should have received a
//  copy of the license along with this program; if not, license is  
//  available online at <http://qlnet.sourceforge.net/License.html>.
//   
//  QLNet is a based on QuantLib, a free-software/open-source library
//  for financial quantitative analysts and developers - http://quantlib.org/
//  The QuantLib license is available online at http://quantlib.org/license.shtml.
//  
//  This program is distributed in the hope that it will be useful, but WITHOUT
//  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
//  FOR A PARTICULAR PURPOSE.  See the license for more details.
using System;
using System.Collections.Generic;

namespace QLNet
{
   //! calibration helper for Equity models
   public class VanillaOptionHelper : CalibrationHelper
   {
      public VanillaOptionHelper( Period maturity,
                                Calendar calendar,
                                double s0,
                                double strikePrice,
                                Handle<Quote> volatility,
                                Handle<YieldTermStructure> riskFreeRate,
                                Handle<YieldTermStructure> dividendYield,
                                CalibrationHelper.CalibrationErrorType errorType = CalibrationErrorType.RelativePriceError)
         :base(volatility, riskFreeRate, errorType)
      {
         maturity_ = maturity;
         calendar_ = calendar;
         s0_ = new Handle<Quote>(new SimpleQuote(s0));
         strikePrice_ = strikePrice;
         dividendYield_ = dividendYield;

         dividendYield.registerWith(update);
      }

      public VanillaOptionHelper( Period maturity,
                                Calendar calendar,
                                Handle<Quote> s0,
                                double strikePrice,
                                Handle<Quote> volatility,
                                Handle<YieldTermStructure> riskFreeRate,
                                Handle<YieldTermStructure> dividendYield,
                                CalibrationHelper.CalibrationErrorType errorType = CalibrationErrorType.RelativePriceError)
         :base(volatility, riskFreeRate, errorType)
      {
         maturity_ = maturity;
         calendar_ = calendar;
         s0_ = s0;
         strikePrice_ = strikePrice;
         dividendYield_ = dividendYield;
          
         s0.registerWith(update);
         dividendYield.registerWith(update);
      }

      public override void addTimesTo(List<double> t) {}

      protected override void performCalculations()
      {
         exerciseDate_ = calendar_.advance(termStructure_.link.referenceDate(), maturity_);
         tau_ = termStructure_.link.timeFromReference(exerciseDate_);
         type_ = strikePrice_ * termStructure_.link.discount(tau_) >=
                        s0_.link.value() * dividendYield_.link.discount(tau_)
                    ? Option.Type.Call
                    : Option.Type.Put;
        StrikedTypePayoff payoff = new PlainVanillaPayoff(type_, strikePrice_);
        Exercise exercise = new EuropeanExercise(exerciseDate_);
        option_ = new VanillaOption(payoff, exercise);
        base.performCalculations();
      }

      public override double modelValue()
      {
         calculate();
         option_.setPricingEngine(engine_);
         return option_.NPV();
      }

      public override double blackPrice(double volatility)
      {
         calculate();
         double stdDev = volatility * Math.Sqrt(maturity());
         return Utils.blackFormula( type_, strikePrice_ * termStructure_.link.discount(tau_),
            s0_.link.value() * dividendYield_.link.discount(tau_), stdDev);
      }

      public double maturity()  { calculate(); return tau_; }
      
   
      private Period maturity_;
      private Calendar calendar_;
      private Handle<Quote> s0_;
      private double strikePrice_;
      private Handle<YieldTermStructure> dividendYield_;
      private Date exerciseDate_;
      private double tau_;
      private Option.Type type_;
      private VanillaOption option_;
      
      public class Export
      {
         int k = 0;
         bool headerPop_ = false;
         System.IO.FileStream file_;
         System.IO.StreamWriter writer_;
         VanillaOptionHelper option_;
         public Export(string path, VanillaOptionHelper option) :
            this(path)
         {
            option_ = option;
         }
         public Export(string path)
         {
            file_ = new System.IO.FileStream(path, System.IO.FileMode.Create);
            writer_ = new System.IO.StreamWriter(file_);
         }
         public string Header()
         {
            headerPop_ = true;
            return "strike\tmaturity\tvolatility";
         }
         public void DoExport()
         {
            if (!headerPop_) writer_.WriteLine(Header());
            writer_.WriteLine(this.ToString());
            writer_.Flush();
         }
         public void Iterate()
         {
            file_.Seek(0, System.IO.SeekOrigin.Begin);
            k++;
         }
         public void setItem(VanillaOptionHelper option) { option_ = option; }
         public override string ToString()
         {
            return option_.strikePrice_ + "\t" + option_.maturity_.length() + "\t" + option_.volatility_.link.value();
         }
      }
    }

   
}
