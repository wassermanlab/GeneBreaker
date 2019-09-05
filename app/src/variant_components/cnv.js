import React from 'react';

function CNV(props) {
    //start, end, copy_change
    if (props.type !== "cnv") {
      return null;
    }
    return (
      <React.Fragment>
        <label>Start position</label>
        <input type="number" className="form-control" min="1"
          name={"var" + props.var + "_start"}
          value={props.start}
          onChange={props.handleInputChange} />
        <label>End position</label>
        <input type="number" className="form-control" min="2"
          name={"var" + props.var + "_end"}
          value={props.end}
          onChange={props.handleInputChange} />
        <label>Copy Change</label>
        <input type="number" className="form-control" min="-1"
          name={"var" + props.var + "_copy_change"}
          value={props.copy_change}
          onChange={props.handleInputChange} />
        <small className="form-text text-muted">
          input an integer copy number change where -1 is a deletion and positive numbers greater than 1 indicated multiplications.
        </small>
      </React.Fragment>
    )
  }


export default CNV;