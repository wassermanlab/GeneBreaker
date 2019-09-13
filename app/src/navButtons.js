import React from 'react';

function NavButtons(props) {
  if (props.page === 1) {
    return (<div><button type="button" className="btn btn-primary float-right" onClick={props.next}>Variant 1</button></div>)
  }
  else if (props.page === 2) {
    return (
      <div>
        <button type="button" className="btn btn-primary float-right" onClick={props.get_vars}>Design Family</button>
        <button type="button" className="btn btn-primary float-right" onClick={props.next}>Variant 2</button>
        <button type="button" className="btn btn-primary float-right" onClick={props.back}>Variant 1</button> </div>)
  }
  else if (props.page === 3) {
    return (
      <div>
        <button type="button" className="btn btn-primary float-right" onClick={props.get_vars}>Design Family</button>
        <button type="button" className="btn btn-primary float-right" onClick={props.back}>Variant 1</button> </div>)
  }
  return null;
}

export default NavButtons;
