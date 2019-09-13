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

  function option(item, index){
    if (props.type === "transcript") {
      return <option key={item.qualifiers.uid} value={item.qualifiers.uid} chrom={item.chrom}> {item.qualifiers.name} </option>
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