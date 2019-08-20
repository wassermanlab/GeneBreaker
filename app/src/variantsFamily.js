import React from 'react';
import { makeStyles } from '@material-ui/core/styles';
import Button from '@material-ui/core/Button';
import Typography from '@material-ui/core/Typography';
import Grid from '@material-ui/core/Grid';

const useStyles = makeStyles(theme => ({
  button: {
    margin: theme.spacing(1),
    float: "right"
  },
  input: {
    display: 'none',
  },
}));

export default function Family(props) {
  // props: 
  //   number
  //   genome
  //   sex
  //   gene_uid 

  const classes = useStyles();

  return (
    <React.Fragment>
      <Grid item xs={12}>
      <Typography variant="h5" component="h2">
              Family 
            </Typography>
            </Grid>
      <Grid item xs={12}>
      <Button variant="contained" color="primary" className={classes.button}>
        Get VCF
      </Button>
      <Button variant="contained" className={classes.button}>
        Back
      </Button>
      </Grid>
      </React.Fragment>
    );
}