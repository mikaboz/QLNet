/*
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
using System.Collections.Generic;
using System.Linq;

namespace QLNet
{
   /// <summary>
   /// Interface pour faciliter la manipulation des modèles à deux facteurs (par rapport aux multi-factors en général)
   /// </summary>
   /*
   public static class TwoFactorModelExtension
   {
      public static double[,] GenerateCorrelations(double? rho)
      {
         if (rho == null) return null;
         else return new double[,]
         {
            {1,(double)rho },
            {(double)rho,1 }
         };
      }
      public static IEnumerable<T> GeneratFactors<T>(T first, T second)
         where T : OneFactorModel
      {
         return new T[] { first, second }.AsEnumerable();

      }
      public static T1 First<T1,T2>(this ITwoFactorModel<T1,T2> model)
         where T1:OneFactorModel
         where T2:OneFactorModel
      {
         return (T1)((MultiFactorModel)model).factors_[0];
      }
      public static T2 Second<T1,T2>(this ITwoFactorModel<T1,T2> model)
         where T1 : OneFactorModel
         where T2 : OneFactorModel
      {
         return (T2)((MultiFactorModel)model).factors_[1];
      }
      public static Parameter Rho<T1,T2>(this ITwoFactorModel<T1,T2> model)
         where T1:OneFactorModel
         where T2 : OneFactorModel
      {
         return ((MultiFactorModel)model).Cor(0, 1);
      }
   }
   public interface ITwoFactorModel<T1, T2>
      where T1: OneFactorModel
      where T2: OneFactorModel
   {}
   */

   
   public abstract class TwoFactorModel : ShortRateModel
   {
      #region Attributes
      public bool IsCorrelated { get; set; }
      protected OneFactorModel first_;
      protected OneFactorModel second_;
      #endregion

      #region Constructors
      private static double[,] GenerateCorrelations(double? rho)
         {
            if (rho == null) return null;
            else return new double[,]
            {
            {1,(double)rho },
            {(double)rho,1 }
            };
         }
      private static IEnumerable<T> GeneratFactors<T>(T first, T second)
            where T : OneFactorModel
         {
            return new T[] { first, second }.AsEnumerable();
      }
      protected TwoFactorModel(OneFactorModel first, OneFactorModel second, double? rho = null) :
         base(GeneratFactors(first, second))
      {
         first_ = first;
         second_ = second;
         int index = 0;
         for (int k = 0; k < first_.Arguments.Count; k++)
            arguments_[index++] = first_.Arguments[k];
         for (int k = 0; k < second_.Arguments.Count; k++)
            arguments_[index++] = second_.Arguments[k];

         if (rho == null)
            IsCorrelated = false;
         else
         {
            arguments_.Add(new ConstantParameter((double)rho, new BoundaryConstraint(-1.0, 1.0)));
            IsCorrelated = true;
         }
      }
      #endregion

      #region Accessors
      public virtual OneFactorModel First { get { return first_; } }
      public virtual OneFactorModel Second { get { return second_; } }
      protected Parameter rho_
      {
         get
         {
            if (IsCorrelated)
            {
               int initialIndex = first_.Arguments.Count + second_.Arguments.Count - 1; // On somme sur chaque facteur le nombre d'arguments individuels
               return arguments_[initialIndex + 1];
            }
            else return new FixedParameter(0.0);
         }
      }
      protected double Rho { get { return rho_.value(0.0); } }
      #endregion

      #region Dynamics
      public override ShortRateModel.Dynamics dynamics()
      {
         return new TwoFactorModel.Dynamics(this);
      }
      public new class Dynamics : ShortRateModel.Dynamics
      {
         public StochasticProcess1D First { get; private set; }
         public StochasticProcess1D Second { get; private set; }
         public double Rho { get; private set; }

         public Dynamics(TwoFactorModel model)
         {
            First = ((OneFactorModel.Dynamics)model.first_.dynamics()).Process;
            First = ((OneFactorModel.Dynamics)model.second_.dynamics()).Process;
            Rho = model.Rho;
         }
         public override double shortRate(double t, Vector variables)
         {
            return ShortRate(t, variables[0], variables[1]);
         }
         public virtual double ShortRate(double t, double x, double y)
         {
            return x + y;
         }
         public override double variable(double t, double r)
         {
            throw new NotImplementedException("On ne récupère pas la varible quand il y a deux facteurs");
         }
      }
      #endregion

      #region Tree
      //! Recombining two-dimensional tree discretizing the state variable
      public override Lattice tree(TimeGrid grid)
      {
         TwoFactorModel.Dynamics dyn = (TwoFactorModel.Dynamics)dynamics();
         TrinomialTree tree1 = new TrinomialTree(dyn.First, grid);
         TrinomialTree tree2 = new TrinomialTree(dyn.Second, grid);
         return (Lattice)(new ShortRateTree(tree1, tree2, dyn));
      }
      public class ShortRateTree : TreeLattice2D<ShortRateTree, TrinomialTree>, IGenericLattice
      {
         protected override ShortRateTree impl()
         {
            return this;
         }
         Dynamics dynamics_;
         //! Plain tree build-up from short-rate dynamics
         public ShortRateTree(TrinomialTree tree1, TrinomialTree tree2, TwoFactorModel.Dynamics dynamics)
            : base(tree1, tree2, dynamics.Rho)
         {
            dynamics_ = dynamics;
         }
         public double discount(int i, int index)
         {
            int modulo = tree1_.size(i);
            int index1 = index % modulo;
            int index2 = index / modulo;

            double x = tree1_.underlying(i, index1);
            double y = tree2_.underlying(i, index2);

            double r = dynamics_.ShortRate(timeGrid()[i], x, y);
            return Math.Exp(-r * timeGrid().dt(i));
         }
         #region Interface
         public double underlying(int i, int index)
         {
            throw new NotImplementedException();
         }
         #endregion
      }
      #endregion
   }

   public abstract class TwoFactorAffineModel : TwoFactorModel, IAffineShortRateModel
   {
      #region Constructors
      protected TwoFactorAffineModel(OneFactorModel first, OneFactorModel second, double? rho = null) :
         base(first, second, rho)
      { }
      #endregion

      #region IAffineModel implémentation
      public abstract double A(double t, double T);
      public Vector Bvect(double t, double T)
      {
         Vector v = new Vector(2)
         {
            [0] = ((OneFactorAffineModel)first_).B(t, T),
            [1] = ((OneFactorAffineModel)second_).B(t, T)
         };
         return v;
      }
      public abstract double DiscountBondOption(Option.Type type, double strike, double maturity, double bondMaturity);
      public virtual double Discount(double t)
      {
         OneFactorModel.Dynamics Dyn1 = (OneFactorModel.Dynamics)First.dynamics();
         OneFactorModel.Dynamics Dyn2 = (OneFactorModel.Dynamics)Second.dynamics();
         Vector vec = new Vector(2)
         {
            [0] = Dyn1.shortRate(0.0, Dyn1.Process.x0()),
            [1] = Dyn2.shortRate(0.0, Dyn2.Process.x0())
         };
         return DiscountBond(0.0, t, vec);
      }
      public virtual double DiscountBond(double t, double T, Vector factors)
      {
         return A(t, T) * Math.Exp(-(Bvect(t, T) * factors));
      }
      public double DiscountBond(double t, double T, double xRate, double yRate)
      {
         Vector vec = new Vector(2)
         {
            [0] = xRate,
            [1] = yRate
         };
         return DiscountBond(t, T, vec);
      }
      #endregion



   }
}

