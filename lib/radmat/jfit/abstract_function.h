#ifndef ABSTRACT_FUNCTION_H
#define ABSTRACT_FUNCTION_H

#include <vector>
#include <string>

namespace jFit
{

  template<typename RType,      // return type
    typename PType,             // parameter type <vector>
    unsigned short NP,          // npars
    typename IType,             // input type <vector>
    unsigned short NI >         // ninputs
      struct AbstractFitFunctionNames
      {
        AbstractFitFunction(void) {}
        ~AbstractFitFunction(void) {};

        virtual RType operator()(const std::vector<IType> &) const = 0;
        virtual std::string getFitType(void) const = 0; 

        unsigned short getNPars(void) const {return NP;}
        unsigned short getNUnfixedPars(void) const; 



        // PARAMETER NAMING
        void setParNames(const std::vector<std::string> &names);
        void setParName(const int parNum, const std::string &name);
        std::vector<std::string> getParNames() const { return parNames; };
        std::string getParName(const int parNum) const;
        int getParNum(const std::string &name) const;

        // DEFAULT VALUES (e.g. fit start values)
        void setDefaultParValues(const std::vector<PType> &values);
        void setDefaultParValue(const int parNum, const PType &value);
        void setDefaultParValue(const std::string &name, const double &value);
        const std::vector<> getDefaultParValues() const { return defaultParValues; };
        double getDefaultParValue(int parNum) const;
        double getDefaultParValue(string name) const;

        // DEFAULT ERRORS (e.g. fit start errors)
        void setDefaultParErrors(vector<double> errors);
        void setDefaultParError(int parNum, double value);
        void setDefaultParError(string name, double value);
        vector<double> getDefaultParErrors() const { return defaultParErrors; };
        double getDefaultParError(int parNum) const;
        double getDefaultParError(string name) const;

        // PARAM FIXING
        void fixParam(int parNum);
        void fixParam(string name);
        void releaseParam(int parNum);
        void releaseParam(string name);

        bool isParamFixed(int parNum) const;
        bool isParamFixed(string name) const;

        // PARAM RANGE FIXING
        void setParamUpperLimit(int parNum, double value);
        void setParamUpperLimit(string name, double value);
        void setParamLowerLimit(int parNum, double value);
        void setParamLowerLimit(string name, double value);
        void setParamLimits(int parNum, double low, double high);
        void setParamLimits(string name, double low, double high);

        void releaseParamLimits(int parNum);
        void releaseParamLimits(string name);

        double getParamUpperLimit(int parNum) const;
        double getParamUpperLimit(string name) const;
        double getParamLowerLimit(int parNum) const;
        double getParamLowerLimit(string name) const;

        bool isParamUpperLimited(int parNum) const;
        bool isParamUpperLimited(string name) const;
        bool isParamLowerLimited(int parNum) const;
        bool isParamLowerLimited(string name) const;

        private:
        int nPars;
        vector<string> parNames;
        vector<double> defaultParValues;
        vector<double> defaultParErrors;
        vector<bool> parFixed;

        vector<double> parLowLimits;
        vector<double> parHighLimits;
        vector<bool> parUpperLimited;
        vector<bool> parLowerLimited;


      };



} // jFit

#endif /* ABSTRACT_FUNCTION_H */
