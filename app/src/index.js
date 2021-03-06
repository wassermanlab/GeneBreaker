import React from 'react';
import ReactDOM from 'react-dom';
import 'bootstrap/dist/css/bootstrap.css';
import 'bootstrap/dist/js/bootstrap.js';
import AppRouter from './app'
import * as serviceWorker from './serviceWorker';
import './index.css'
  

// ReactDOM.render(
//     <VForm />,
//     document.getElementById('root')
//   );

ReactDOM.render(
  <AppRouter />,
  document.getElementById('root')
);

// If you want your app to work offline and load faster, you can change
// unregister() to register() below. Note this comes with some pitfalls.
// Learn more about service workers: https://bit.ly/CRA-PWA
serviceWorker.unregister();
