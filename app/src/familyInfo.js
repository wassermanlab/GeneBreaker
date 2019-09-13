import React from 'react';

class FInfo extends React.Component {
  constructor(props) {
    super(props);
    this.state = {
    };
  }

  render() {
    if (this.props.page !== 4) {
      return null;
    }
    let row2;
    if (this.props.vars.var1.proband === "1/1") {
      row2 = (<td colSpan="4">SAME AS V1</td>)
    }
    else if (this.props.vars.var2 === "") {
      row2 = (<td colSpan="4">NONE</td>)
    }
    else {
      row2 = [<td>{this.props.vars.var2.chrom}</td>,
      <td>{this.props.vars.var2.pos}</td>,
      <td>{this.props.vars.var2.ref}</td>,
      <td>{this.props.vars.var2.alt}</td>]
    }
    return (
      <React.Fragment>
      {/* variants */}
      <h2>Variants</h2>
      <table className="table">
        <thead>
          <tr>
            <th scope="col">*</th>
            <th scope="col">CHROM</th>
            <th scope="col">POS</th>
            <th scope="col">REF</th>
            <th scope="col">ALT</th>
          </tr>
        </thead>
        <tbody>
          <tr>
            <th scope="row">V1</th>
            <td>{this.props.vars.var1.chrom}</td>
            <td>{this.props.vars.var1.pos}</td>
            <td>{this.props.vars.var1.ref}</td>
            <td>{this.props.vars.var1.alt}</td>
          </tr>
          <tr>
            <th scope="row">V2</th>
            {row2}
          </tr>
        </tbody>
      </table>
     {/* family */}
      <h2>Family</h2>
      </React.Fragment>
    )
  }
}

export default FInfo;
