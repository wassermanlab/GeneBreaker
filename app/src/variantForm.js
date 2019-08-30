import React from 'react';
import { Formik} from 'formik';
import PageGeneral from './page_general';
import PageVariant from './page_variant';

const sleep = ms => new Promise(resolve => setTimeout(resolve, ms));

const required = value => (value ? undefined : 'Required');

class Wizard extends React.Component {
  static Page = ({ children }) => children;

  constructor(props) {
    super(props);
    this.state = {
      page: 0,
      values: props.initialValues,
    };
  }

  next = values =>
    this.setState(state => ({
      page: Math.min(state.page + 1, this.props.children.length - 1),
      values,
    }));

  previous = () =>
    this.setState(state => ({
      page: Math.max(state.page - 1, 0),
    }));

  validate = values => {
    const activePage = React.Children.toArray(this.props.children)[
      this.state.page
    ];
    return activePage.props.validate ? activePage.props.validate(values) : {};
  };

  handleSubmit = (values, bag) => {
    const { children, onSubmit } = this.props;
    const { page } = this.state;
    const isLastPage = page === React.Children.count(children) - 1;
    if (isLastPage) {
      return onSubmit(values, bag);
    } else {
      bag.setTouched({});
      bag.setSubmitting(false);
      this.next(values);
    }
  };

  render() {
    const { children } = this.props;
    const { page, values } = this.state;
    const activePage = React.Children.toArray(children)[page];
    const isLastPage = page === React.Children.count(children) - 1;
    return (
      <Formik
        initialValues={values}
        enableReinitialize={false}
        validate={this.validate}
        onSubmit={this.handleSubmit}
        render={({ values, handleSubmit, isSubmitting, handleReset }) => (
          <form onSubmit={handleSubmit}>
            {activePage}
            <div className="buttons">
              {page > 0 && (
                <button
                  type="button"
                  className="secondary"
                  onClick={this.previous}
                >
                  « Previous
                </button>
              )}

              {!isLastPage && <button type="submit">Next »</button>}
              {isLastPage && (
                <button type="submit" disabled={isSubmitting}>
                  Submit
                </button>
              )}
            </div>

          </form>
        )}
      />
    );
  }
}


const VFrom = () => (
  <div className="VFrom">
    <h1>Variant Simulator</h1>
    <Wizard
      initialValues={{
        gene_uid: "",
        genome: "hg38",
        chr: "chrY",
        sex: "",
        var1: {
          type: "",
          region: "",
          zygosity: "",
          impact: {
            clinvar: {
              id: ""
            },
            cnv: {
              start: "",
              end: "",
              copy_change: ""
            },
            indel: {
              length: "",
              start: ""
            },
            mei: {
              element: "",
              start: ""
            },
            snv: {
              snv_type: "",
              start: ""
            },
            str: {
              start: "",
              end: "",
              motif: "",
              length: ""
            }
          }
        },
        var2: {
          type: "",
          region: "",
          zygosity: "",
          impact: {
            clinvar: {
              id: ""
            },
            cnv: {
              start: "",
              end: "",
              copy_change: ""
            },
            indel: {
              length: "",
              start: ""
            },
            mei: {
              element: "",
              start: ""
            },
            snv: {
              snv_type: "",
              start: ""
            },
            str: {
              start: "",
              end: "",
              motif: "",
              length: ""
            }
          }
        }

      }}
      onSubmit={(values, actions) => {
        sleep(300).then(() => {
          window.alert(JSON.stringify(values, null, 2));
          actions.setSubmitting(false);
        });
      }}
    >
      <Wizard.Page>
        <PageGeneral required={required} />
      </Wizard.Page>
      <Wizard.Page>
        <PageVariant var={"var1"} required={required} />
      </Wizard.Page>
      <Wizard.Page>
        <PageVariant var={"var2"} required={required} />
      </Wizard.Page>
    </Wizard>
  </div>
);

export default VFrom;
