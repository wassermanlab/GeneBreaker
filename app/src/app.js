import React from "react";
import { BrowserRouter as Router, Route } from "react-router-dom";
import VForm from './variantForm'

function Index() {
  return <h2>Home</h2>;
}

function AppRouter() {
  return (
    <Router>
      <div>
        <Route path="/" exact component={Index} />
        <Route path="/variants/" component={VForm} />
        {/* <Route path="/family/" component={Family} /> */}
      </div>
    </Router>
  );
}

export default AppRouter;
