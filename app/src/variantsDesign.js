import React from 'react';
import Row from 'react-bootstrap/Row';
import Col from 'react-bootstrap/Col';
import Container from 'react-bootstrap/Container';
import Form from 'react-bootstrap/Form';
import General from './variantsGeneral'
import Family from './variantsFamily'
import Variant from './variantsVariant'

class DesignVariants extends React.Component {

  constructor(props) {
    super(props);
    this.state = {
      currentStep: 0,         // panel options: 0=generalInfo, 1=var1, 2=var2, 3=family
      genome: "hg38",
      sex: "XX",
      chromosome: "",
      gene_uid: "",
      var1: {
        type: "",
        region: "",
        zygosity: "",
        impact: ""
      },
      var2: ""
    };
    // Bind the submission to handleChange() 
    this.handleChange = this.handleChange.bind(this)
  }

  // Use the submitted data to set the state
  handleChange(event) {
    const { name, value } = event.target
    this.setState({
      [name]: value
    })
}

handleSubmit = (event) => {
  event.preventDefault()
  const { email, username, password } = this.state
  alert(`Your registration detail: \n 
    Email: ${email} \n 
    Username: ${username} \n
    Password: ${password}`)
}


render() {
  return (
    <React.Fragment>
      <Container>
        <Form onSubmit={this.handleSubmit}>
          {/* general info */}
          <General currentStep={this.state.currentStep} />
          {/* variant 1 */}
          <Variant currentStep={this.state.currentStep} />
          {/* variant 2 */}
          <Variant currentStep={this.state.currentStep} />
          {/* family info */}
          <Family currentStep={this.state.currentStep} />

        </Form>
      </Container>
    </React.Fragment>
  );
}
}

export default DesignVariants; 