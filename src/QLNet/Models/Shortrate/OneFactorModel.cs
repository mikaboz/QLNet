﻿/*
 Copyright (C) 2008 Siarhei Novik (snovik@gmail.com)
 Copyright (C) 2010 Philippe Real (ph_real@hotmail.com)

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
using System.Linq;

namespace QLNet
{
   //! Single-factor short-rate model abstract class
   /*! \ingroup shortrate */

   public abstract class OneFactorModel : ShortRateModel
   {
      protected OneFactorModel(int nArguments) : base(nArguments)
      { }

      //! Base class describing the short-rate dynamics
      public new class Dynamics : ShortRateModel.Dynamics
      {
         public StochasticProcess1D Process { get; set; }
         public Dynamics(StochasticProcess1D process)
         {
            Process = process;
         }

         //! Compute state variable from short rate
         //public abstract double variable(double t, double rate);
         public override double shortRate(double t, Vector variables)
         {
            return shortRate(t, variables[0]);
         }
         //! Compute short rate from state variable
         public virtual double shortRate(double t, double variable)
         {
            return variable;
         }
         public override double variable(double t, double r)
         {
            return r;
         }
      }

      //! Return by default a trinomial recombining tree
      public override Lattice tree(TimeGrid grid)
      {
         OneFactorModel.Dynamics dyn = (OneFactorModel.Dynamics)dynamics();
         TrinomialTree trinomial = new TrinomialTree(dyn.Process, grid,false);
         return new ShortRateTree(trinomial, dyn, grid);
      }

      //! Recombining trinomial tree discretizing the state variable
      public class ShortRateTree : TreeLattice1D<ShortRateTree>, IGenericLattice
      {
         protected override ShortRateTree impl()
         {
            return this;
         }

         //! Plain tree build-up from short-rate dynamics
         public ShortRateTree(TrinomialTree tree, Dynamics dynamics, TimeGrid timeGrid) :
            base(timeGrid, tree.size(1))
         {
            tree_ = tree;
            dynamics_ = dynamics;
         }
         //! Tree build-up + numerical fitting to term-structure
         public ShortRateTree(TrinomialTree tree, Dynamics dynamics,TermStructureFittingParameter.NumericalImpl theta, TimeGrid timeGrid) :
            base(timeGrid, tree.size(1))
         {
            tree_ = tree;
            dynamics_ = dynamics;
            theta.reset();
            double value = 1.0;
            double vMin = -100.0;
            double vMax = 100.0;
            for (int i = 0; i < (timeGrid.size() - 1); i++)
            {
               double discountBond = theta.termStructure().link.discount(t_[i + 1]);
               Helper finder = new Helper(i, discountBond, theta, this);
               Brent s1d = new Brent();
               s1d.setMaxEvaluations(1000);
               value = s1d.solve(finder, 1e-7, value, vMin, vMax);
               theta.change(value);
            }
         }

         public int size(int i)
         {
            return tree_.size(i);
         }

         public double discount(int i, int index)
         {
            double x = tree_.underlying(i, index);
            double r = dynamics_.shortRate(timeGrid()[i], x);
            return Math.Exp(-r * timeGrid().dt(i));
         }

         public override double underlying(int i, int index)
         {
            return tree_.underlying(i, index);
         }

         public int descendant(int i, int index, int branch)
         {
            return tree_.descendant(i, index, branch);
         }

         public double probability(int i, int index, int branch)
         {
            return tree_.probability(i, index, branch);
         }

         private TrinomialTree tree_;
         private Dynamics dynamics_;

         public class Helper : ISolver1d
         {
            private int size_;
            private int i_;
            private Vector statePrices_;
            private double discountBondPrice_;
            private TermStructureFittingParameter.NumericalImpl theta_;
            private ShortRateTree tree_;

            public Helper(int i,
               double discountBondPrice,
               TermStructureFittingParameter.NumericalImpl theta,
               ShortRateTree tree)
            {
               size_ = tree.size(i);
               i_ = i;
               statePrices_ = tree.statePrices(i);
               discountBondPrice_ = discountBondPrice;
               theta_ = theta;
               tree_ = tree;
               theta_.setvalue(tree.timeGrid()[i], 0.0);
            }

            public override double value(double theta)
            {
               double value = discountBondPrice_;
               theta_.change(theta);
               for (int j = 0; j < size_; j++)
                  value -= statePrices_[j] * tree_.discount(i_, j);
               return value;
            }
         }
      }
   }

   public abstract class OneFactorAffineModel : OneFactorModel, IAffineShortRateModel
   {
      #region Constructor
      protected OneFactorAffineModel(int nArguments)
         : base(nArguments)
      {}
      #endregion
      #region IAffineModel implémentation et adaptation 1D
      public abstract double A(double t, double T);
      public Vector Bvect(double t, double T)
      {
         double[] b = new double[] { B(t, T) };
         return new Vector(b);
      }
      public abstract double B(double t, double T);
      public virtual double DiscountBondOption(Option.Type type, double strike, double maturity, double bondMaturity)
      {
         throw new NotImplementedException();
      }
      public double Discount(double t)
      {
         OneFactorModel.Dynamics Dynamics = (OneFactorModel.Dynamics)dynamics();
         Vector vec = new Vector(1)
         {
            [0] = Dynamics.shortRate(0.0, Dynamics.Process.x0())
         };
         return ((IAffineShortRateModel)this).DiscountBond(0.0, t, vec);
      }
      public double DiscountBond(double t, double T, Vector factors)
      {
         return A(t, T) * Math.Exp(-(Bvect(t, T) * factors));
      }
      public double DiscountBond(double t, double T, double rate)
      {
         Vector v = new Vector(1)
         {
            [0] = rate
         };
         return DiscountBond(t, T, v);
      }
      #endregion
   }
}
