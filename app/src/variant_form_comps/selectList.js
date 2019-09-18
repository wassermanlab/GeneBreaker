import React, { useState, useEffect } from 'react';

function SelectList(props) {
  // PROPS: 
  //// title -> title of your list 
  //// url -> url to make request to 
  //// name -> name 
  //// value -> value
  //// type -> for option
  //// watch -> something to watch
  //// onChange -> change handler 

  const [list, setList] = useState([]);

  useEffect(() => {

    async function fetchUrl() {
      const response = await fetch(props.url);
      const json = await response.json();
      setList(json);
    }
    fetchUrl();
  }, [props.url]);

  function option(item, index) {
    if (props.type === "transcript") {
      return <option key={item.qualifiers.uid} value={item.qualifiers.uid} chrom={item.chrom}> {item.qualifiers.name} </option>
    } else if (props.type === "clingen") {
      return (<option key={index + item.qualifiers.uid} value={item.qualifiers.uid} start={item.start + 1} end={item.end} copy_change={item.qualifiers.copy_number_change}>
        Start: {item.start + 1} &emsp;&emsp; End: {item.end} &emsp;&emsp; Clinical assertion: {item.qualifiers.clinical_assertion} &emsp;&emsp; Copy change: {item.qualifiers.copy_number_change} </option>
      )
    } else if (props.type === 'str') {
      return (<option key={index + item.qualifiers.uid} value={item.qualifiers.uid}>
        Start: {item.start + 1} &emsp;&emsp; End: {item.end} &emsp;&emsp; Motif: {item.qualifiers.motif} &emsp;&emsp; Pathogenicity: {item.qualifiers.pathogenicity} </option>)
    } else if (props.type === 'clinvar') {
      return (<option key={index + item.qualifiers.uid} value={item.qualifiers.uid}>
        Start: {item.start + 1} &emsp;&emsp; Ref: {item.qualifiers.ref} &emsp;&emsp; Alt: {item.qualifiers.alt} &emsp;&emsp; Clinical Significance: {item.qualifiers.CLNSIG} </option>
      )
    }
  }

  return (
    <React.Fragment>
      <div className="form-group">
        <label>{props.title}</label>
        <select className="form-control"
          name={props.name}
          value={props.value}
          onChange={props.onChange}
          size="5">
          <option key="0" value=""></option>
          {list.map((item, index) => (
            option(item, index)
          ))}
        </select>
      </div>
    </React.Fragment>
  )
}

export default SelectList;