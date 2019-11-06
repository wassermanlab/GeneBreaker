import React from "react";
import { BrowserRouter as Router, Route } from "react-router-dom";
import MasterForm from "./v2/masterForm";
import Home from "./home";

function AppRouter() {
  return (
    <Router>
        <Route path="/" exact component={Home} />
        <Route path="/variants/" component={MasterForm} />
    </Router>
  );
}

export default AppRouter;
