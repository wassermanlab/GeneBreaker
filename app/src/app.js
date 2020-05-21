import React from "react";
import { BrowserRouter as Router, Route } from "react-router-dom";
import MasterForm from "./v2/masterForm";
import Home from "./home";
import MoreInfo from "./more_info";
import PremadeCases from "./premade_cases";

function AppRouter() {
  return (
    <Router>
        <Route path="/" exact component={Home} />
        <Route path="/variants/" component={MasterForm} />
        <Route path="/premade_cases/" component={PremadeCases} />
        <Route path="/more_info/" component={MoreInfo} />
    </Router>
  );
}

export default AppRouter;
