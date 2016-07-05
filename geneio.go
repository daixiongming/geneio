// Package geneio provides interfaces for gene I/O functions.
package geneio

import (
	"github.com/biogo/biogo/feat/gene"

	"io"
)

// Reader is the common gene.Interface reader interface.
type Reader interface {
	// Read reads a gene.Interface, returning any error that occurs during the
	// read.
	Read() (gene.Interface, error)
}

// Writer is the common gene.Interface writer interface.
type Writer interface {
	// Write writes a gene.Interface, returning the number of bytes written
	// and any error that occurs during the write.
	Write(gene.Interface) (n int, err error)
}

// Scanner wraps a Reader to provide a convenient loop interface for reading
// feature data.  Successive calls to the Scan method will step through the
// genes of the provided Reader. Scanning stops unrecoverably at EOF or the
// first error.
type Scanner struct {
	r   Reader
	g   gene.Interface
	err error
}

// NewScanner returns a Scanner to read from r.
func NewScanner(r Reader) *Scanner {
	return &Scanner{r: r}
}

// NewScannerFromFunc returns a Scanner to read genes returned by calls to f.
func NewScannerFromFunc(f func() (gene.Interface, error)) *Scanner {
	return NewScanner(funcReader(f))
}

// Next advances the Scanner past the next gene, which will then be
// available through the Gene method. It returns false when the scan stops,
// either by reaching the end of the input or an error. After Next returns
// false, the Error method will return any error that occurred during
// scanning, except that if it was io.EOF, Error will return nil.
func (s *Scanner) Next() bool {
	if s.err != nil {
		return false
	}
	s.g, s.err = s.r.Read()
	return s.err == nil
}

// Error returns the first non-EOF error that was encountered by the Scanner.
func (s *Scanner) Error() error {
	if s.err == io.EOF {
		return nil
	}
	return s.err
}

// Gene returns the most recent gene read by a call to Next.
func (s *Scanner) Gene() gene.Interface {
	return s.g
}

// funcReader is a wrapper type that implements Reader.
type funcReader func() (gene.Interface, error)

// Read calls f and returns any gene.Interface and error returned by f.
func (f funcReader) Read() (gene.Interface, error) {
	return f()
}
