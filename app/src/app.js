import React from "react";
import { BrowserRouter as Router, Route } from "react-router-dom";
import MasterForm from "./variant_form_comps/masterForm";
import MasterForm2 from "./v2/masterForm";
import Home from "./home";

function AppRouter() {
  return (
    <Router>
        <Route path="/" exact component={Home} />
        <Route path="/variants/" component={MasterForm} />
        <Route path="/v2/" component={MasterForm2} />
        {/* <Route path="/family/" component={Family} /> */}
    </Router>
  );
}

export default AppRouter;
