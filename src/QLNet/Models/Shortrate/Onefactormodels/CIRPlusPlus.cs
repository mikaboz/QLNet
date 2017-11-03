/*
 Copyright (C) 2017 Mikael BOZON (mikael.bozon@gmail.com)
  
 This file is part of QLNet Project https://github.com/amaggiulli/qlnet

 QLNet is free software: you can redistribute it and/or modify it
 under the terms of the QLNet license.  You should have received a
 copy of the license along with this program; if not, license is  
 available online at <http://qlnet.sourceforge.net/License.html>.
  
 QLNet is a based on QuantLib, a free-software/open-source library
 for financial quantitative analysts and developers - http://quantlib.org/
 The QuantLib license is available online at http://quantlib.org/license.shtml.
 
 This program is distributed in the hope that it will be useful, but WITHOUT
 ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 FOR A PARTICULAR PURPOSE.  See the license for more details.
*/

using System;

namespace QLNet
{
   public class CIRPlusPlus : CoxIngersollRoss, ITermStructureConsistentModel
   {
      public double ValueForFitting(double t)
      {
         double h = Math.Sqrt(Kappa * Kappa + 2 * Sigma * Sigma);
         double expth = Math.Exp(t * h);
         double temp = (Kappa + h) * (expth - 1);
         double numerator1 = 2 * Kappa * Theta * (expth - 1);
         double denominator1 = 2 * h + temp;
         double numerator2 = 4 * h * h * expth;
         double denominator2_ = 2 * h + temp;
         double denominator2 = denominator2_ * denominator2_;
         return numerator1 / denominator1 + r0_ * numerator2 / denominator2;
      }
      public class FittingParameter : TermStructureFittingParameter
      {
         public FittingParameter(CIRPlusPlus model) :
            base(new Impl(model))
         { }
         private new class Impl : TermStructureFittingParameter.Impl
         {
            public Impl(CIRPlusPlus model) :
               base(model)
            {}
            public override double value(Vector p, double t)
            {
               return model_.TermStructureInitialForwardRate(t) - ((CIRPlusPlus)model_).ValueForFitting(t);
            }
         }
      }
      private FittingParameter phi_;

      public CIRPlusPlus(Handle<YieldTermStructure> termStructure, double r0, double kappa = 0.1, double theta = 0.1, double sigma = 0.1) :
         base(r0,kappa,theta,sigma)
      {
         termStructure_ = termStructure;
         termStructure_.registerWith(update);
         generateArguments();
      }

      protected override void generateArguments()
      {
         phi_ = new FittingParameter(this);
      }
      #region ITermStructureConsistentModel
      private Handle<YieldTermStructure> termStructure_;
      public Handle<YieldTermStructure> TermStructure { get { return termStructure_; }}
      #endregion
      
      public override double A(double t, double T)
      {
         double P0T = termStructure_.link.discount(T);
         double P0t = termStructure_.link.discount(t);
         double A0T = base.A(0, T);
         double A0t = base.A(0, t);
         double AtT = base.A(t, T);
         double B0T = base.B(0, T);
         double B0t = base.B(0, t);
         double BtT = base.B(t, T);
         double numerator = P0T * A0t * Math.Exp(-B0t * r0_);
         double denominator = P0t * A0T * Math.Exp(-B0T * r0_);
         return (numerator / denominator) * AtT * Math.Exp(BtT) * phi_.value(t);
      }
      /// <summary>
      /// Incertain
      /// </summary>
      public override double DiscountBondOption(Option.Type type, double strike, double maturity, double bondMaturity)
      {
         // double t = 0 (=t dans les formules d'évaluation) par rapport à la date de référence.
         double P0S = termStructure_.link.discount(bondMaturity);
         double P0T = termStructure_.link.discount(maturity);
         double P0t = termStructure_.link.discount(0.0);
         double A0S = base.A(0, bondMaturity);
         double A0T = base.A(0, maturity);
         double A0t = base.A(0, 0.0);
         double B0S = base.B(0, bondMaturity);
         double B0T = base.B(0, maturity);
         double B0t = base.B(0, 0.0);
         double psiCIRStrikeNumerator = P0T * A0S * Math.Exp(-B0S * r0_);
         double psiCIRStrikeDenominator = P0S * A0T * Math.Exp(-B0T * r0_);
         double psiCIRStrike = psiCIRStrikeNumerator / psiCIRStrikeDenominator;

         // qu'est-ce que x dans psiCIR(t,T,tau,X,x,alpa)?? Cf Brigo Mercurio p.103
         return base.DiscountBondOption(type, psiCIRStrike, maturity, bondMaturity);
      }
      public override ShortRateModel.Dynamics dynamics()
      {
         return new CIRPlusPlus.Dynamics(phi_,new SquareRootProcess(r0_,Kappa,Theta,Sigma));
      }
      public new class Dynamics : OneFactorModel.Dynamics
      {
         private FittingParameter fitting_;
         public Dynamics(FittingParameter fitting, SquareRootProcess process) :
            base(process)
         {
            fitting_ = fitting;
         }
         public override double variable(double t, double r)
         {
            return r - fitting_.value(t);
         }
         public override double shortRate(double t, double x)
         {
            return x + fitting_.value(t);
         }
      }
   }
}
