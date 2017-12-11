using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;


using QLNet;
namespace PDD
{
   public class TESTTT { }
   public class MCBarriereEngine<RNG, S> : McSimulation<SingleVariate, RNG, S>, IGenericEngine
      where RNG : IRSG, new()
      where S : IGeneralStatistics, new()
   {
      protected GeneralizedBlackScholesProcess process_;
      protected int? timeSteps_, timeStepsPerYear_;
      protected int? requiredSamples_, maxSamples_;
      protected double? requiredTolerance_;
      protected bool brownianBridge_;
      protected ulong seed_;

      public MCBarriereEngine(GeneralizedBlackScholesProcess process,
                                 int? timeSteps,
                                 int? timeStepsPerYear,
                                 bool brownianBridge,
                                 bool antitheticVariate,
                                 bool controlVariate,
                                 int? requiredSamples,
                                 double? requiredTolerance,
                                 int? maxSamples,
                                 ulong seed)
         : base(antitheticVariate, controlVariate)
      {
         process_ = process;
         timeSteps_ = timeSteps;
         timeStepsPerYear_ = timeStepsPerYear;
         requiredSamples_ = requiredSamples;
         maxSamples_ = maxSamples;
         requiredTolerance_ = requiredTolerance;
         brownianBridge_ = brownianBridge;
         seed_ = seed;

         Utils.QL_REQUIRE(timeSteps != null || timeStepsPerYear != null, () => "no time steps provided");
         Utils.QL_REQUIRE(timeSteps == null || timeStepsPerYear == null, () =>
                   "both time steps and time steps per year were provided");
         if (timeSteps != null)
            Utils.QL_REQUIRE(timeSteps > 0, () => "timeSteps must be positive, " + timeSteps + " not allowed");
         if (timeStepsPerYear != null)
            Utils.QL_REQUIRE(timeStepsPerYear > 0, () =>
           "timeStepsPerYear must be positive, " + timeStepsPerYear + " not allowed");

         process_.registerWith(update);
      }

      public virtual void calculate()
      {
         base.calculate(requiredTolerance_, requiredSamples_, maxSamples_);
         results_.value = mcModel_.sampleAccumulator().mean();
         if (FastActivator<RNG>.Create().allowsErrorEstimate != 0)
            results_.errorEstimate = mcModel_.sampleAccumulator().errorEstimate();
      }

      // McSimulation implementation
      protected override TimeGrid timeGrid()
      {
         Date lastExerciseDate = arguments_.exercise.lastDate();
         double t = process_.time(lastExerciseDate);
         if (timeSteps_ != null)
         {
            return new TimeGrid(t, timeSteps_.Value);
         }
         else if (timeStepsPerYear_ != null)
         {
            int steps = (int)(timeStepsPerYear_ * t);
            return new TimeGrid(t, Math.Max(steps, 1));
         }
         else
         {
            Utils.QL_FAIL("time steps not specified");
            return null;
         }
      }

      protected override IPathGenerator<IRNG> pathGenerator()
      {
         int dimensions = process_.factors();
         TimeGrid grid = timeGrid();
         IRNG generator = (IRNG)FastActivator<RNG>.Create().make_sequence_generator(dimensions * (grid.size() - 1), seed_);
            return new PathGenerator<IRNG>(process_, grid, generator, brownianBridge_);
      }

      protected override double? controlVariateValue()
      {
         throw new NotImplementedException();
      }

      protected override PathPricer<IPath> pathPricer()
      {
         PlainVanillaPayoff payoff = arguments_.payoff as PlainVanillaPayoff;
         if (payoff == null)
            throw new NotSupportedException("Plain vanilla payoff needed");

         return new BarrierPathPricer(process_.blackVolatility(),
            payoff.optionType(),
            payoff.strike(),
            arguments_.barrierType,
            arguments_.barrier.Value,
            arguments_.rebate.Value,
            process_.riskFreeRate().link.discount(arguments_.exercise.lastDate()));
      }

      #region PricingEngine
      protected BarrierOption.Arguments arguments_ = new BarrierOption.Arguments();
      protected BarrierOption.Results results_ = new BarrierOption.Results();

      public IPricingEngineArguments getArguments() { return arguments_; }
      public IPricingEngineResults getResults() { return results_; }
      public void reset() { results_.reset(); }

      #region Observer & Observable
      // observable interface
      private readonly WeakEventSource eventSource = new WeakEventSource();
      public event Callback notifyObserversEvent
      {
         add { eventSource.Subscribe(value); }
         remove { eventSource.Unsubscribe(value); }
      }

      public void registerWith(Callback handler) { notifyObserversEvent += handler; }
      public void unregisterWith(Callback handler) { notifyObserversEvent -= handler; }
      protected void notifyObservers()
      {
         eventSource.Raise();
      }

      public void update() { notifyObservers(); }
      #endregion
      #endregion
   }


   public class BarrierPathPricer : PathPricer<IPath>
   {
      PlainVanillaPayoff payOff;
      Barrier.Type barrierType_;
      double barrier_;
      double rebate_;
      double discount_;
      Handle<BlackVolTermStructure> sigma_;
      public BarrierPathPricer(Handle<BlackVolTermStructure> sigma, Option.Type optionType, double strike, Barrier.Type barrierType, double barrier, double rebate, double discount)
      {
         sigma_ = sigma;
         payOff = new PlainVanillaPayoff(optionType, strike);
         barrierType_ = barrierType;
         barrier_ = barrier;
         rebate_ = rebate;
         discount_ = discount;
      }

      private bool triggered(double forwardPrice)
      {
         switch (barrierType_)
         {
            case Barrier.Type.DownIn:
            case Barrier.Type.DownOut:
               return forwardPrice < barrier_;
            case Barrier.Type.UpIn:
            case Barrier.Type.UpOut:
               return forwardPrice > barrier_;
            default:
               Utils.QL_FAIL("unknown type");
               return false;
         }
      }

      public double value(IPath p)
      {
         Path path = (Path)p;
         Vector values = path.values();
         if (triggered(values.Max()))
            return rebate_;
         else
            if (new BrownianBridgeExtremumCondition(barrier_ ,barrierType_,sigma_).MaxReached(path))
            return rebate_;
         else
            return payOff.value(path.back()) * discount_;

      }

      public class BrownianBridgeExtremumCondition
      {
         MersenneTwisterUniformRng rng = new MersenneTwisterUniformRng();
         private enum Extremum  {Minimum, Maximum}
         private Extremum extremum_;
         private Handle<BlackVolTermStructure> sigma_;
         private double barrier_;
         public BrownianBridgeExtremumCondition(double barrier, Barrier.Type barrierType, Handle<BlackVolTermStructure> sigma)
         {
            barrier_ = barrier;
            sigma_ = sigma;
            switch (barrierType)
            {
            case Barrier.Type.DownIn:
            case Barrier.Type.DownOut:
                  extremum_ = Extremum.Minimum;
                  break;
            case Barrier.Type.UpIn:
            case Barrier.Type.UpOut:
                  extremum_ = Extremum.Maximum;
                  break;
            default:
                  throw new NotSupportedException("unknown type");
            }
         }
         public bool MaxReached(Path path)
         {
            Vector values = path.values();
            bool maxReached = false;
            TimeGrid grid = path.timeGrid();
            MersenneTwisterUniformRng URNG = new MersenneTwisterUniformRng();
            int i = 0;
            while(i<values.Count-1 && maxReached == false)
            {
               double t = grid[i];
               double T = grid[i + 1];
               double h = 1 / (double)grid.dt(i);
               double xt = values[i];
               double xT = values[i + 1];
               double diff = xT - xt;
               double sigma2 = sigma_.link.blackForwardVariance(t, T, 0)/(T-t);
               double unif = URNG.nextReal();
               switch (extremum_)
               {
                  case Extremum.Maximum:
                     double pk = Math.Exp(2 * (barrier_ - xt) * (barrier_ - xT) / (h * sigma2));
                     if (unif >= pk)
                        maxReached = true;
                     break;
                  case Extremum.Minimum:
                     throw new NotImplementedException("not yet implemented");
                  default:
                     throw new NotSupportedException("Extremum type not known");
               }
               i++;
            }
            return maxReached;

         }
      }


   }

}
