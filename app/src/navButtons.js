import React from 'react';

function NavButtons(props) {
  if (props.page === 1) {
    return (<div><button type="button" className="btn btn-primary float-right" onClick={props.next}>Next</button></div>)
  }
  else if (props.page === 2) {
    return (
      <div>
        <button type="button" className="btn btn-primary float-right" onClick={props.next}>Next</button>
        <button type="button" className="btn btn-primary float-right" onClick={props.back}>Back</button> </div>)
  }
  else if (props.page === 3) {
    return (
      <div>
        <button type="button" className="btn btn-primary float-right" onClick={props.back}>Back</button> </div>)
  }
  return null;
}

export default NavButtons;
