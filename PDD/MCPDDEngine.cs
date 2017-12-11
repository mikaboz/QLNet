using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

using QLNet;
namespace PDD
{
   public class MCPDDEngine<RNG, S> : McSimulation<SingleVariate, RNG, S>, IGenericEngine
            where RNG : IRSG, new()
            where S : IGeneralStatistics, new()
   {
      // data members
      protected GeneralizedBlackScholesProcess process_;
      private int stepsPerConditionalPeriod_;
      private int? requiredSamples_;
      private int? maxSamples_;
      private double? requiredTolerance_;
      private bool brownianBridge_;
      ulong seed_;

      bool startWithConditonalPeriod_ = false;

      // constructor
      public MCPDDEngine(
           GeneralizedBlackScholesProcess process,
           int stepsPerConditionalPeriod,
           bool brownianBridge,
           bool antitheticVariate,
           bool controlVariate,
           int? requiredSamples,
           double? requiredTolerance,
           int? maxSamples,
           ulong seed) : base(antitheticVariate, controlVariate)
      {
         process_ = process;
         stepsPerConditionalPeriod_ = stepsPerConditionalPeriod;
         requiredSamples_ = requiredSamples;
         maxSamples_ = maxSamples;
         requiredTolerance_ = requiredTolerance;
         brownianBridge_ = brownianBridge;
         seed_ = seed;
         process_.registerWith(update);
      }

      public void calculate()
      {
         base.calculate(requiredTolerance_, requiredSamples_, maxSamples_);
         results_.value = this.mcModel_.sampleAccumulator().mean();
         if (FastActivator<RNG>.Create().allowsErrorEstimate != 0)
            results_.errorEstimate =
                this.mcModel_.sampleAccumulator().errorEstimate();
      }


      

      protected override TimeGrid timeGrid()
      {
         Calendar calendar = process_.blackVolatility().link.calendar();
         DayCounter dc = process_.blackVolatility().link.dayCounter();
         Date evaluationDate = calendar.adjust(Settings.evaluationDate());
         Period conditionalPeriod = arguments_.ConditionalPeriod;
         Period provPeriod = arguments_.ProvFrequency;
         Date maturityDate = calendar.adjust(arguments_.exercise.lastDate());
         
         List<double> times = new List<double>();
         times.Add(dc.yearFraction(evaluationDate, maturityDate));
         Date upperDate = maturityDate;
         Date lowerDate = calendar.advance(upperDate, -conditionalPeriod);
         if (lowerDate <= evaluationDate)
         {
            double conditionalPeriodTime = dc.yearFraction(evaluationDate, maturityDate);
            double dt = conditionalPeriodTime / (double)stepsPerConditionalPeriod_;
            for (int i = 0; i < stepsPerConditionalPeriod_; i++)
               times.Add(times.Last() - dt);
            return new TimeGrid(times);//PDDTimeGrid(times, stepsPerConditionalPeriod_, true);
         }
         while (lowerDate > evaluationDate)
         {
            double conditionalPeriodTime = dc.yearFraction(lowerDate, upperDate);
            double dt = conditionalPeriodTime / (double)stepsPerConditionalPeriod_;
            for (int i = 0; i < stepsPerConditionalPeriod_; i++)
               times.Add(times.Last() - dt);
            upperDate = lowerDate;
            lowerDate = calendar.advance(upperDate, conditionalPeriod - provPeriod);
            if (lowerDate <= evaluationDate)
               break;
            double NoConditionalPeriodTime = dc.yearFraction(lowerDate, upperDate);
            times.Add(times.Last() - NoConditionalPeriodTime);
         }
         times.Add(0);
         times.Reverse();
         startWithConditonalPeriod_ = true;
         return new TimeGrid(times);
      }
      /*
      protected override TimeGrid timeGrid()
      {
         Console.WriteLine("Timegrid building");
         Calendar calendar = process_.blackVolatility().link.calendar();
         DayCounter dc = process_.blackVolatility().link.dayCounter();
         Date evaluationDate = calendar.adjust(Settings.evaluationDate());
         Period conditionalPeriod = arguments_.ConditionalPeriod;
         Period provPeriod = arguments_.ProvFrequency;
         Date maturityDate = calendar.adjust(arguments_.exercise.lastDate());

         List<TimeGrid> timeGrids = new List<TimeGrid>();

         Date upperDate = maturityDate;
         Date lowerDate = calendar.advance(upperDate, -conditionalPeriod);
         if (lowerDate < evaluationDate)
            lowerDate = evaluationDate;
         double conditionalPeriodTime = dc.yearFraction(lowerDate, upperDate);
         timeGrids.Add(new TimeGrid(conditionalPeriodTime, stepsPerConditionalPeriod_));
         if (lowerDate == evaluationDate)
            return timeGrids.First();

         // Il va faloir construire une sequence
         double firstPeriodTime = 0;
         List<double> phasis = new List<double>();
         while (lowerDate > evaluationDate)
         {
            upperDate = lowerDate;
            lowerDate = calendar.advance(upperDate, conditionalPeriod - provPeriod);
            Console.ReadLine();
            if (lowerDate < evaluationDate)
            {
               lowerDate = evaluationDate;
               firstPeriodTime = dc.yearFraction(lowerDate, upperDate);
               break;
            }
            phasis.Add(dc.yearFraction(lowerDate, upperDate));
            upperDate = lowerDate;
            lowerDate = calendar.advance(upperDate, -conditionalPeriod);
            if (lowerDate < evaluationDate)
               lowerDate = evaluationDate;
            conditionalPeriodTime = dc.yearFraction(lowerDate, upperDate);
            timeGrids.Add(new TimeGrid(conditionalPeriodTime, stepsPerConditionalPeriod_));
         }

         return new SequenceTimeGrid(timeGrids, firstPeriodTime, phasis);
      }
      */
      /*
      protected override TimeGrid timeGrid()
      {
         Console.WriteLine("Timegrid building");
         Calendar calendar = process_.blackVolatility().link.calendar();
         DayCounter dc = process_.blackVolatility().link.dayCounter();
         Date evaluationDate = calendar.adjust(Settings.evaluationDate());
         Period conditionalPeriod = arguments_.ConditionalPeriod;
         Period provPeriod = arguments_.ProvFrequency;
         Date maturityDate = calendar.adjust(arguments_.exercise.lastDate());

         List<TimeGrid> timeGrids = new List<TimeGrid>();

         Date upperDate = maturityDate;
         Date lowerDate = calendar.advance(upperDate, -conditionalPeriod);
         if (lowerDate < evaluationDate)
            lowerDate = evaluationDate;
         double conditionalPeriodTime = dc.yearFraction(lowerDate, upperDate);
         timeGrids.Add(new TimeGrid(conditionalPeriodTime, stepsPerConditionalPeriod_));
         if (lowerDate == evaluationDate)
            return timeGrids.First();

         // Il va faloir construire une sequence
         double firstPeriodTime = 0;
         List<double> phasis = new List<double>();
         while (lowerDate>evaluationDate)
         {
            upperDate = lowerDate;
            lowerDate = calendar.advance(upperDate, conditionalPeriod - provPeriod);
            Console.ReadLine();
            if (lowerDate < evaluationDate)
            {
               lowerDate = evaluationDate;
               firstPeriodTime = dc.yearFraction(lowerDate, upperDate);
               break;
            }
            phasis.Add(dc.yearFraction(lowerDate,upperDate));
            upperDate = lowerDate;
            lowerDate = calendar.advance(upperDate, -conditionalPeriod);
            if (lowerDate < evaluationDate)
               lowerDate = evaluationDate;
            conditionalPeriodTime = dc.yearFraction(lowerDate, upperDate);
            timeGrids.Add(new TimeGrid(conditionalPeriodTime, stepsPerConditionalPeriod_));
         }
         return new SequenceTimeGrid(timeGrids, firstPeriodTime, phasis);
      }

      */
      // McSimulation implementation
      /*
      protected override TimeGrid timeGrid()
      {
         List<double> times_ = new List<double>();

         Calendar calendar = process_.blackVolatility().link.calendar();
         Date evaluationDate = calendar.adjust(Settings.evaluationDate());
         Period conditionalPeriod = arguments_.ConditionalPeriod;
         Period provPeriod = arguments_.ProvFrequency;
         Date maturityDate = calendar.adjust(arguments_.exercise.lastDate());

         bool firstDateIsConditionalPeriod = false;
         List<Date> dates = new List<Date>();
         dates.Add(maturityDate);
         do
         {
            Date date = calendar.advance(dates.Last(), -conditionalPeriod);
            if (date < evaluationDate)
            {
               firstDateIsConditionalPeriod = true;
               dates.Add(evaluationDate);
            }
            else
            {
               dates.Add(date);
               Date date2 = calendar.advance(dates.Last(), provPeriod-conditionalPeriod);
               if (date < evaluationDate)
                  dates.Add(evaluationDate);
               else
                  dates.Add(date2);
            }
         }
         while (dates.Last() > evaluationDate);
         dates.Reverse();

         DayCounter dc = process_.blackVolatility().link.dayCounter();
         bool tempIsConditionalPeriod = firstDateIsConditionalPeriod;

         times_.Add(0.0);
         for (int i = 0; i < dates.Count-1; i++)
         {
            double periodTime = dc.yearFraction(dates[i], dates[i+1]);
            if (tempIsConditionalPeriod)
            {
               double dt = periodTime / (double)stepsPerConditionalPeriod_;
               for (int k = 0; k < stepsPerConditionalPeriod_; k++)
                  times_.Add(times_.Last() + dt);
            }
            else
            {
               times_.Add(times_.Last() + periodTime);
            }
            tempIsConditionalPeriod = !tempIsConditionalPeriod;
         }
         return new TimeGrid(times_, times_.Count);
      }
      */

      protected override IPathGenerator<IRNG> pathGenerator()
      {
         TimeGrid grid = this.timeGrid();
         IRNG gen = (IRNG)new RNG().make_sequence_generator(grid.size() - 1, seed_);
         return new PathGenerator<IRNG>(process_, grid,
                                        gen, brownianBridge_);
      }

      protected override double? controlVariateValue()
      {
         IPricingEngine controlPE = this.controlPricingEngine();
         Utils.QL_REQUIRE(controlPE != null, () => "engine does not provide control variation pricing engine");

         PDD.Arguments controlArguments =
                 (PDD.Arguments)controlPE.getArguments();
         controlArguments = arguments_;
         controlPE.calculate();

         PDD.Results controlResults =
             (PDD.Results)(controlPE.getResults());

         return controlResults.value;

      }

      protected override PathPricer<IPath> pathPricer()
      {
         PDD.Payoff payoff = arguments_.payoff as PDD.Payoff;
         Utils.QL_REQUIRE(payoff != null, () => "non-PDD payoff given");

         GeneralizedBlackScholesProcess process = process_ as GeneralizedBlackScholesProcess;
         Utils.QL_REQUIRE(process != null, () => "Black-Scholes process required");

         return new PDDPathPricer(payoff, process.riskFreeRate().link.discount(timeGrid().Last()),stepsPerConditionalPeriod_,startWithConditonalPeriod_);
      
      }

      #region PricingEngine
      protected PDD.Arguments arguments_ = new PDD.Arguments();
      protected PDD.Results results_ = new PDD.Results();

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
}
