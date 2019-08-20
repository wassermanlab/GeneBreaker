import React from 'react';
import Paper from '@material-ui/core/Paper';
import Grid from '@material-ui/core/Grid';
import Typography from '@material-ui/core/Typography';
import CssBaseline from '@material-ui/core/CssBaseline';
import Variant from './variantsVariant'
import General from './variantsGeneral'
import Family from './variantsFamily'
import Container from '@material-ui/core/Container';

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
  }
  render() {
    return (
      <React.Fragment>
        <CssBaseline />
        <Container maxWidth="md">
          <Paper>
            <div  style={{ padding: 20 }}>
              <Grid container spacing={3}>
                <Grid item xs={12}>
                  <Typography variant="h4" component="h1" align='center'>
                    Design a Case
            </Typography>
                </Grid>
                {/* stepper */}
                {/* step 0: general information */}
                <General sex="XX" />
                {/* step 1: var 1 */}
                <Variant number="1" genome="hg38" sex="XX" gene_uid="1234" />
                {/* step 2: var 2 */}
                <Variant number="2" genome="hg38" sex="XX" gene_uid="1234" />
                {/* step 3: family */}
                <Family sex="XX" />
              </Grid>
            </div>
          </Paper>
        </Container>
      </React.Fragment>
    );
  }
}

export default DesignVariants; 