import React from 'react';

function Errors(props) {
  return (
    <React.Fragment>{
      props.errors.map((item, index) => (
        <div key={"error_" + index} className="alert alert-danger" role="alert"> {item} </div>))}
    </React.Fragment>
  )
}
export default Errors;

