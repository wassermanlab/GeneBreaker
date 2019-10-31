import React from 'react';
import InputComp from './inputComp'
import SelectComp from './selectComp'
import SelectList from './selectListComp'
import { host } from '../host'

function FamilyInfo(props) {
  if (props.page !== 4) {
    return null;
  }

  return (
    <React.Fragment>
      <h1>Family Info</h1>
     
    </React.Fragment >
  )
}


export default FamilyInfo;