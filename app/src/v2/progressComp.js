import React from 'react';

function Progress(props) {
  let text;
  let value;
  switch (props.page) {
    case 1:
      text = "General Info";
      value = "25";
      break
    case 2:
      text = "Variant 1";
      value = "50";
      break
    case 3:
      text = "Variant 2";
      value = "75";
      break
    case 4:
      text = "Family Info";
      value = "100";
      break
    default:
      break
  }
  const style = value + "%"
  const progress = (
  <React.Fragment>
  <div className="progress" style={{height: "30px", marginBottom: "2%"}}>
    <div className="progress-bar" role="progressbar" style={{width: style, fontSize: "large"}} aria-valuenow={value} aria-valuemin="0" aria-valuemax="100">{text}</div>
  </div>
  </React.Fragment>
  )

  return progress;
}
export default Progress;
