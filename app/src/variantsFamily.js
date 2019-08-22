import React from 'react';
import Row from 'react-bootstrap/Row';
import Col from 'react-bootstrap/Col';
import Button from 'react-bootstrap/Button';

class Family extends React.Component {
  render() {
    if (this.props.currentStep !== 4) { // Prop: The current step
      return null
    }
    // The markup for the Step 4 UI
    return(
      <div>
        <h2>Family Information</h2>
        <Row>
          <Col sm={10}></Col>
          <Col sm={1}>
            <Button variant="secondary">Back</Button>
          </Col>
          <Col sm={1}>
            <Button variant="primary">Submit</Button>
          </Col>
        </Row>
      </div>
    )
  }
}


export default Family; 