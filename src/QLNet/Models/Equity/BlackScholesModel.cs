using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace QLNet
{
   public interface IParametricVolatility
   {
      double IntegratedSquareValue(double t, double T);
   }
   public class ConstantVolatilityParameter : ConstantParameter, IParametricVolatility
   {
      public ConstantVolatilityParameter(double constantVol) :
         base(constantVol, new PositiveConstraint())
      { }
      public double IntegratedSquareValue(double t, double T)
      {
         double vol = value(0);
         return vol * vol * (T-t);
      }
   }
   public abstract class TimeDeterministParameter : Parameter, IParametricVolatility
   {
      private GlobalPositivityConstraint globalPositivityConstraint_;
      public Period horizon_;
      public void SetHorizon(double horizon)
      {
         if (globalPositivityConstraint_ != null)
            globalPositivityConstraint_.SetHorizon(horizon);
      }
      public double? Horizon
      {
         get
         {
            return globalPositivityConstraint_?.Horizon;
         }
      }
      protected TimeDeterministParameter(Impl impl, Period horizon, bool needGlobalPositivity)
      {
         impl_ = impl;
         params_ = impl.params_;
         horizon_ = horizon;
         
         if (needGlobalPositivity)
         {
            globalPositivityConstraint_ = new GlobalPositivityConstraint((Impl)impl_, horizon);
            constraint_ = globalPositivityConstraint_;
         }
         else constraint_ = new NoConstraint();
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
         public override double value(Vector p, double t)
         {
            CheckAndAdaptValue(ref t, p);
            return valueImpl(p, t);
         }
         public virtual void CheckAndAdaptValue(ref double t, Vector p) {  }
         public abstract double valueImpl(Vector p, double t);
         public virtual void CheckAndAdaptIntegrability(ref double t, ref double T, Vector p) { }
         public double IntegratedSquareValue(Vector p, double t, double T)
         {
            CheckAndAdaptIntegrability(ref t, ref T, p);
            return IntegratedSquareValueImpl(p, t, T);
         }
         public abstract double IntegratedSquareValueImpl(Vector p, double t, double T);
      }
      private class GlobalPositivityConstraint : Constraint
      {
         public double? Horizon { get { return ((Impl)impl_).horizon_; } }
         public void SetHorizon(double horizon)
         {
            ((Impl)impl_).SetHorizon(horizon);
         }
         public GlobalPositivityConstraint(TimeDeterministParameter.Impl functionImpl, Period Horizon) :
            base(new Impl(functionImpl))
         { }
         public class Impl : IConstraint
         {
            public double? horizon_;
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
               if (horizon_ == null)
                  throw new ArgumentNullException("null horizon, must be set from termStructure and more specificly from ParametricVolatilityTermStructure");
               double horizon = (double)horizon_;
               Vector initialValue = new Vector(1, horizon / 2);
               Constraint boundaryConstraint = new BoundaryConstraint(0, horizon);
               Problem problem = new Problem(new CostFunctionImpl(functionImpl_, parameter), boundaryConstraint, initialValue);
               EndCriteria endCriteria = new EndCriteria(100, null, 10e-4, 10e-4, null );
               LevenbergMarquardt lm = new LevenbergMarquardt();
               lm.minimize(problem, endCriteria);
               double minimum = problem.value(problem.currentValue());
               bool test = minimum > 0;
               Console.WriteLine("GlobalPositivityConstraint: " + test);
               return test;
            }
           
            public class CostFunctionImpl : CostFunction
            {
               public TimeDeterministParameter.Impl function_;
               public Vector parameters_;
               public CostFunctionImpl(TimeDeterministParameter.Impl function, Vector parameters)
               {
                  function_ = function;
                  parameters_ = parameters;
               }
               public override double value(Vector x)
               {
                  double t = x[0];
                  double value = function_.value(parameters_,t);
                  return value;
               }
               public override Vector values(Vector x)
               {
                  Vector v = new Vector(1)
                  {
                     [0] = value(x)
                  };
                  return v;
               }
            }
         }
      }
   }
   /// <summary>
   /// if c belongs to [-inf,-1[ : polynomial
   /// if c = 1 : linear
   /// if c belongs to ]-1,0[ : sqrt
   /// if = 0 : constr
   /// if c belongs to ]0,1[ : inverse and 0 - integrable
   /// ifc belongs to [1,+inf[ : inverse and not-0-integrable : TO PREVENT
   /// </summary>
   public class InverseVolatilityParameter : TimeDeterministParameter
   {
      public enum FormConstraint { NoSpecific, Polynomial, Inverse, Sqrt }
      public InverseVolatilityParameter(double a, double b, double c, Period horizon, FormConstraint formConstraint = FormConstraint.NoSpecific) :
         base(new Impl(a,b,c),horizon,true)
      {
         Constraint integrability = new IntegrabilityConstraintImpl(this.params_);
         Constraint form = new FormConstraintImpl(formConstraint, params_);
         Constraint composite = new CompositeConstraint(integrability, form);
         constraint_ = new CompositeConstraint(constraint_, composite);
      }
      public bool CheckParameters(FormConstraint formConstraint)
      {
         double c = params_[2];
         switch (formConstraint)
         {
            case FormConstraint.Inverse:
               return c > 0 & c < 1;
            case FormConstraint.Polynomial:
               return c <= 1;
            case FormConstraint.Sqrt:
               return c >= 1 && c <= 0;
            case FormConstraint.NoSpecific:
               return true;
            default:
               throw new NotSupportedException("not support form Constraint in checking parameters");
         }
      }
      
      public new class Impl : TimeDeterministParameter.Impl
      {
         public Impl(double a, double b, double c) :
            base(3)
         {
            params_[0] = a;
            params_[1] = b;
            params_[2] = c;
         }

         public override double IntegratedSquareValueImpl(Vector p, double t, double T)
         {
            double A = p[0];
            double B = p[1];
            double C = p[2];
            if (C >= 1)
               throw new NotSupportedException("c into inverseVolatilityParameter doesn't belong to ]-inf,1[ : integrability constraint did not impact");
            //double t = (t_ == 0) ? (t_ + 0.0001) : t_;
            //double T = (T_ == 0) ? (T_ + 0.0001) : T_;
            return A * (T - t) + 2 * B * A * (Math.Pow(T, 1 - C) - Math.Pow(t, 1 - C)) / (1 - C) + B * B * (Math.Pow(T, 1 - 2 * C) - Math.Pow(t, 1 - 2 * C)) / (1 - 2 * C);
         }
         public override void CheckAndAdaptValue(ref double t, Vector p)
         {
            double c = p[2];
            t = (t == 0&& c>=0) ? t + QLNet.Const.QL_EPSILON : t;
         }
         public override double valueImpl(Vector p, double t)
         {
            double A = p[0];
            double B = p[1];
            double C = p[2];
            return A + B / Math.Pow(t, C);
         }
      }
      public class FormConstraintImpl : ProjectedIndividualConstraint
      {
         private static Constraint SelectedConstraint(FormConstraint form)
         {
            switch (form)
            {
               case FormConstraint.Inverse:
                  return new BoundaryConstraint(0, double.MaxValue);
               case FormConstraint.Polynomial:
                  return new BoundaryConstraint(double.MinValue, -1);
               case FormConstraint.Sqrt:
                  return new BoundaryConstraint(-1, 0);
               default:
                  throw new NotSupportedException("Not support formConstraint");
            }
         }
         public FormConstraintImpl(FormConstraint formConstraint, Vector parameters) :
            base(SelectedConstraint(formConstraint),parameters,2)
         { }
         public override bool test(Vector p)
         {
            bool test = base.test(p);
            Console.WriteLine("formConstraint: " + test);
            return test;
         }
      }
      public class IntegrabilityConstraintImpl : ProjectedIndividualConstraint
      {
         public IntegrabilityConstraintImpl(Vector parameters) :
            base(new ExtendedBoundaryConstraint(double.MinValue,1, true, true),parameters,2)
         { }
         public override bool test(Vector p)
         {
            bool test = base.test(p);
            Console.WriteLine("IntegrabilityConstraint: " + test);
            return test;
         }
      }
   }
   public class ExponentialVolatilityParameter : TimeDeterministParameter
   {
      public ExponentialVolatilityParameter(double A, double alpha, double B, double beta, Period horizon) :
         base(new Impl(A,alpha,B,beta),horizon,true)
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

         public override double IntegratedSquareValueImpl(Vector p, double t, double T)
         {
            double A = p[0];
            double alpha = p[1];
            double B = p[2];
            double beta = p[3];
            return A * A * (Math.Exp(-2 * alpha * t) - Math.Exp(-2 * alpha * T)) / (2 * alpha)
                + 2 * A * B * (Math.Exp(-(alpha + beta) * t) - Math.Exp(-(alpha + beta) * T)) / (alpha + beta)
                + B * B * (Math.Exp(-2 * beta * t) - Math.Exp(-2 * beta * T)) / (2 * beta);
         }
         public override double valueImpl(Vector p, double t)
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
      public Period Horizon
      {
         get
         {
            if (parameter_ is TimeDeterministParameter timeDeterministParameter)
               return timeDeterministParameter.horizon_;
            else return null;
         }
      }
      public Parameter parameter_;
      public ParametricVolatilityTermStructure(IParametricVolatility p, Date referenceDate, DayCounter dayCounter, Calendar calendar, BusinessDayConvention bdc = BusinessDayConvention.Following) :
         base(referenceDate,calendar,bdc,dayCounter)
      {
         parameter_ = (Parameter)p;
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
         return ((IParametricVolatility)parameter_).IntegratedSquareValue(0, t);
      }
   }
   
   public class BlackScholesModel : CalibratedModel
   {
      public GeneralizedBlackScholesProcess process_;
      public BlackScholesModel(Handle<Quote> spot, Handle<YieldTermStructure> dividendTS,  Handle<YieldTermStructure> riskFreeTS, ParametricVolatilityTermStructure volatilityTS) :
         base(1)
      {
         arguments_[0] = volatilityTS.parameter_;
         Handle<BlackVolTermStructure> blackVolTS = new Handle<BlackVolTermStructure>(volatilityTS);
         process_ = new GeneralizedBlackScholesProcess(spot, dividendTS, riskFreeTS, blackVolTS);
         blackVolTS.registerWith(update);
      }
   }
}
