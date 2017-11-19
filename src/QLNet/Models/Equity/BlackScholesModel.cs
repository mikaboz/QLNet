using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace QLNet
{
   public class ConstantVolatilityParameter : ConstantParameter
   {
      public ConstantVolatilityParameter(double constantVol) :
         base(constantVol, new PositiveConstraint())
      { }

   }
   public class TimeDeterministParameter : Parameter
   {
      public Period horizon_;
      public double Horizon { get { return ((PositivityConstraint)constraint_).Horizon; } }
      protected TimeDeterministParameter(Impl impl, Period horizon) :
         base(impl.Size, impl, new  PositivityConstraint(impl))
      {
         params_ = impl.params_;
         horizon_ = horizon;
      }
      public double IntegratedSquareValue(double t, double T)
      {
         return ((TimeDeterministParameter.Impl)impl_).IntegratedSquareValue(params_, t, T);
      }
      public abstract new class Impl : Parameter.Impl
      {
         public int Size { get { return params_.Count; } }
         public Vector params_;
         public Impl(int size)
         {
            params_ = new Vector(size);
         }
         public abstract double IntegratedSquareValue(Vector p, double t, double T);
      }
      // test si aux paramètres fixés, pour tout t dans [0,T], sigma(t) > 0
      public class PositivityConstraint : Constraint
      {
         public double Horizon { get { return ((Impl)impl_).horizon_; } }
         public void SetHorizon(double horizon)
         {
            ((Impl)impl_).SetHorizon(horizon);
         }
         public PositivityConstraint(TimeDeterministParameter.Impl functionImpl) :
            base(new Impl(functionImpl))
         { }
         public class Impl : IConstraint
         {
            public double horizon_ = double.MaxValue;
            public TimeDeterministParameter.Impl functionImpl_;
            public Impl(TimeDeterministParameter.Impl functionImpl)
            {
               functionImpl_ = functionImpl;
            }
            public void SetHorizon(double horizon)
            {
               horizon_ = horizon;
            }
            public Vector lowerBound(Vector parameters)
            {
               return new Vector(parameters.size(), Double.MinValue);
            }
            public Vector upperBound(Vector parameters)
            {
               return new Vector(parameters.size(), Double.MaxValue);
            }
            public bool test(Vector parameter)
            {
               Brent brent = new Brent();
               brent.setMaxEvaluations(100);
               //Test si la valeur initiale est négative
               if (functionImpl_.value(parameter, 0) < 0)
               {
                  return false;
               }
               // Alors la volatilité initiale est potivie
               // Dans ce cas, si on trouve un 0 c'est que la volatilité change de signe et passe néfative
               try
               {
                  brent.solve(new Solver(functionImpl_, parameter), 10e-4, horizon_ / 2, 0, horizon_);
               }
               // Si on ne trouve pas de 0, une exception est alors générée
               catch (ArgumentException e)
               {
                  return true;
               }
               return false;
            }
            public class Solver : ISolver1d
            {
               public Vector parameters_;
               public TimeDeterministParameter.Impl functionImpl_;
               public Solver(TimeDeterministParameter.Impl functionImpl, Vector parameters)
               {
                  functionImpl_ = functionImpl;
                  parameters_ = parameters;
               }
               public override double value(double t)
               {
                  return functionImpl_.value(parameters_, t);
               }
            }
         }
      }
      public void SetHorizon(double horizon)
      {
         ((PositivityConstraint)constraint_).SetHorizon(horizon);
      }
   }

   public class InverseVolatilityParameter : TimeDeterministParameter
   {
      public InverseVolatilityParameter(double a, double b, double c, Period horizon) :
         base(new Impl(a,b,c),horizon)
      {}
      public new class Impl : TimeDeterministParameter.Impl
      {
         public Impl(double a, double b, double c) :
            base(3)
         {
            params_[0] = a;
            params_[1] = b;
            params_[2] = c;
         }

         public override double IntegratedSquareValue(Vector p, double t_, double T_)
         {
            double A = p[0];
            double B = p[1];
            double C = p[2];
            double t = (t_ == 0) ? (t_ + 0.0001) : t_;
            double T = (T_ == 0) ? (T_ + 0.0001) : T_;
            return A * (T - t) + 2 * B * A * (Math.Pow(T, 1 - C) - Math.Pow(t, 1 - C)) / (1 - C) + B * B * (Math.Pow(T, 1 - 2 * C) - Math.Pow(t, 1 - 2 * C)) / (1 - 2 * C);
         }
         public override double value(Vector p, double t)
         {
            double t_ = (t == 0) ? (t + 0.0001) : t;
            double A = p[0];
            double B = p[1];
            double C = p[2];
            return A + B / Math.Pow(t_, C);
         }
      }
   }
   public class ExponentialVolatilityParameter : TimeDeterministParameter
   {
      public ExponentialVolatilityParameter(double A, double alpha, double B, double beta, Period horizon) :
         base(new Impl(A,alpha,B,beta),horizon)
      { }
      public new class Impl : TimeDeterministParameter.Impl
      {
         public Impl(double A, double alpha, double B, double beta) :
            base(4)
         {
            params_[0] = A;
            params_[1] = alpha;
            params_[2] = B;
            params_[3] = beta;
         }

         public override double IntegratedSquareValue(Vector p, double t, double T)
         {
            double A = p[0];
            double alpha = p[1];
            double B = p[2];
            double beta = p[3];
            return A * A * (Math.Exp(-2 * alpha * t) - Math.Exp(-2 * alpha * T)) / (2 * alpha)
                + 2 * A * B * (Math.Exp(-(alpha + beta) * t) - Math.Exp(-(alpha + beta) * T)) / (alpha + beta)
                + B * B * (Math.Exp(-2 * beta * t) - Math.Exp(-2 * beta * T)) / (2 * beta);
         }
         public override double value(Vector p, double t)
         {
            double A = p[0];
            double alpha = p[1];
            double B = p[2];
            double beta = p[3];
            return A * Math.Exp(-alpha * t) + B * Math.Exp(-beta * t);
         }
      }
   }

   public class ParametricVolatilityTermStructure : BlackVarianceTermStructure
   {
      public Period Horizon { get
         {
            if (parameter_ is TimeDeterministParameter timeDeterministParameter)
               return timeDeterministParameter.horizon_;
            else return null;
         }
      }
      public Parameter parameter_;
      public ParametricVolatilityTermStructure(Parameter p, Date referenceDate, DayCounter dayCounter, Calendar calendar, BusinessDayConvention bdc = BusinessDayConvention.Following) :
         base(referenceDate,calendar,bdc,dayCounter)
      {
         parameter_ = p;
         if (parameter_ is TimeDeterministParameter timeDetermisitParameter)
            timeDetermisitParameter.SetHorizon(timeFromReference(maxDate()));
      }
      public override Date maxDate()
      {
         if (Horizon is null)
            return Date.maxDate();
         else
            return referenceDate() + Horizon;
      }
      public override double maxStrike()
      {
         return double.MaxValue;
      }
      public override double minStrike()
      {
         return double.MinValue;
      }
      protected override double blackVarianceImpl(double t, double strike)
      {
         if (parameter_ is TimeDeterministParameter timeDeterministParameter)
            return timeDeterministParameter.IntegratedSquareValue(0, t);
         else if (parameter_ is ConstantParameter constantParameter)
         {
            double vol = constantParameter.value(t);
            return vol * vol * t;
         }
         else throw new NotSupportedException("not supported Parameter type");
      }
   }

   public class BlackScholesModel : CalibratedModel
   {
      public BlackScholesMertonProcess process_;
      public BlackScholesModel(Handle<Quote> spot, Handle<YieldTermStructure> dividendTS,  Handle<YieldTermStructure> riskFreeTS, ParametricVolatilityTermStructure volatilityTS) :
         base(1)
      {
         arguments_[0] = volatilityTS.parameter_;
         Handle<BlackVolTermStructure> blackVolTS = new Handle<BlackVolTermStructure>(volatilityTS);
         process_ = new BlackScholesMertonProcess(spot, dividendTS, riskFreeTS, blackVolTS);
         blackVolTS.registerWith(update);
      }
   }
}
