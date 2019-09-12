import React from 'react';

class FForm extends React.Component {
  constructor(props) {
    super(props);
    this.state = {
    };

    this.handleInputChange = this.handleInputChange.bind(this);
  }
 
  handleInputChange(event) {
    const target = event.target;
    const value = target.value;
    const name = target.name;

    // state change for everything else 
    this.setState({
      [name]: value
    },
      () => console.log(this.state));
  }
  render() {
    console.log(this.props.location.state.vcf)
    return (
    <div>HELLO</div>
      )
  }
}

export default FForm;
