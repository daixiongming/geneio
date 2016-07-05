package geneio

import (
	"errors"
	"io"
	"testing"

	"github.com/biogo/biogo/feat/gene"

	"gopkg.in/check.v1"
)

// Hook up gocheck into the "go test" runner.
func Test(t *testing.T) { check.TestingT(t) }

// Create the test suite
type S struct{}

var _ = check.Suite(&S{})

func geneCreator(g gene.Interface, err error) func() (gene.Interface, error) {
	return func() (gene.Interface, error) {
		tmpG, tmpE := g, err
		g, err = nil, io.EOF
		return tmpG, tmpE
	}
}

// Test ReadFromFunc
var readFromFuncTests = []struct {
	f   func() (gene.Interface, error)
	g   gene.Interface
	err string
}{
	{
		geneCreator(&gene.Gene{ID: "A"}, nil),
		&gene.Gene{ID: "A"},
		"",
	},
	{
		geneCreator(nil, errors.New("returning error")),
		nil,
		"returning error",
	},
	{
		geneCreator(&gene.Gene{ID: "B"}, io.EOF),
		&gene.Gene{ID: "B"},
		"",
	},
}

func (s *S) TestReadFromFunc(c *check.C) {
	for _, tt := range readFromFuncTests {
		sc := NewScannerFromFunc(tt.f)

		ok := sc.Next()
		g, err := sc.Gene(), sc.Error()

		// Next() should not return true when there is an error.
		if ok {
			c.Assert(err, check.Equals, nil)
		}

		if tt.err != "" {
			c.Assert(err, check.ErrorMatches, tt.err)
		} else if err != nil {
			c.Errorf("Unexpected error %s", err)
		}

		c.Check(g, check.DeepEquals, tt.g)

		// Second call to Next() should return false.
		c.Assert(sc.Next(), check.Equals, false)
	}
}
