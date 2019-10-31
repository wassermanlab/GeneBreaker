import React from 'react';
// import InputComp from './inputComp'
// import SelectComp from './selectComp'
// import SelectList from './selectListComp'
// import { host } from '../host'

function FamilyInfo(props) {
  if (props.page !== 4) {
    return null;
  }
  console.log(props.vars)
  return (
    <React.Fragment>
      <h1>Family Info</h1>
      <table className="table">
        <thead>
          <tr>
            <th scope="col">Variant</th>
            <th scope="col">Chrom</th>
            <th scope="col">Pos</th>
            <th scope="col">Ref</th>
            <th scope="col">Alt</th>
          </tr>
        </thead>
        <tbody>        
          <tr>
            <th scope="row">1</th>
            <td>test</td>
            <td>test</td>
            <td>test</td>
            <td>test</td>
          </tr>
          <tr>
            <th scope="row">2</th>
            <td>test</td>
            <td>test</td>
            <td>test</td>
            <td>test</td>
          </tr>
        </tbody>
      </table>
    </React.Fragment >
  )
}


export default FamilyInfo;