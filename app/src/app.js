import React from "react";
import { BrowserRouter as Router, Route } from "react-router-dom";
import MasterForm from "./v2/masterForm";
import Home from "./home";
import MoreInfo from "./more_info";
import InheritanceTestingCases from "./inheritance_testing_cases";
import DiseaseScenarios from "./disease_scenarios";
import PatientScenarios from "./patient_scenarios";

function AppRouter() {
  return (
    <Router>
        <Route path="/" exact component={Home} />
        <Route path="/variants/" component={MasterForm} />
        <Route path="/inheritance_testing_cases/" component={InheritanceTestingCases} />
        <Route path="/disease_scenarios/" component={DiseaseScenarios} />
        <Route path="/patient_scenarios/" component={PatientScenarios} />
        <Route path="/more_info/" component={MoreInfo} />
    </Router>
  );
}

export default AppRouter;
