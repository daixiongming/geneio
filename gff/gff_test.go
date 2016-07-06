package gff

import (
	"fmt"
	"io"
	"reflect"
	"strings"
	"testing"

	"github.com/biogo/biogo/feat"
	"github.com/biogo/biogo/io/featio/gff"
	"github.com/go-bio/geneio"
)

// Assert that interfaces are satisfied.
var (
	_ geneio.Reader = (*Reader)(nil)
)

// Test Read
var readTests = []struct {
	Name                   string
	Input                  string
	Error                  string
	GeneTag, TranscriptTag string
	GeneCnt, FeatCnt       int
	IDs, Chrs              []string
	Starts, Ends           []int
	Orientations           []feat.Orientation
}{
	{
		Name: "Normal",
		Input: "" +
			"X\t.\texon\t10\t20\t0\t+\t.\tgene_id A; transcript_id A1;\n" +
			"X\t.\texon\t15\t30\t0\t+\t.\tgene_id A; transcript_id A2;\n" +
			"X\t.\texon\t20\t45\t0\t+\t.\tgene_id A; transcript_id A3;\n" +
			"X\t.\texon\t25\t30\t0\t+\t.\tgene_id A; transcript_id A4;\n" +
			"Y\t.\texon\t50\t90\t0\t-\t.\tgene_id B; transcript_id B1;\n" +
			"Y\t.\tfoo\t100\t200\t0\t-\t.\tgene_id B; transcript_id B1;\n" +
			"Y\t.\tstop_codon\t60\t62\t0\t-\t.\tgene_id B; transcript_id B1;\n" +
			"Y\t.\tstart_codon\t71\t73\t0\t-\t.\tgene_id B; transcript_id B1;\n",
		GeneCnt:      2,
		FeatCnt:      5,
		IDs:          []string{"A", "B"},
		Chrs:         []string{"X", "Y"},
		Starts:       []int{9, 49},
		Ends:         []int{45, 90},
		Orientations: []feat.Orientation{feat.Forward, feat.Reverse},
	},
	{
		Name: "Single forward coding gene",
		Input: "" +
			"X\t.\texon\t2\t70\t0\t+\t.\tgene_id C; transcript_id C1;\n" +
			"X\t.\texon\t80\t90\t0\t+\t.\tgene_id C; transcript_id C1;\n" +
			"X\t.\tstart_codon\t60\t62\t0\t+\t.\tgene_id C; transcript_id C1;\n" +
			"X\t.\tstop_codon\t81\t83\t0\t+\t.\tgene_id C; transcript_id C1;\n",
		GeneCnt:      1,
		FeatCnt:      1,
		IDs:          []string{"C"},
		Chrs:         []string{"X"},
		Starts:       []int{1},
		Ends:         []int{90},
		Orientations: []feat.Orientation{feat.Forward},
	},
	{
		Name: "Single reverse coding gene",
		Input: "" +
			"X\t.\texon\t30\t50\t0\t-\t.\tgene_id D; transcript_id D1;\n" +
			"X\t.\texon\t80\t99\t0\t-\t.\tgene_id D; transcript_id D1;\n" +
			"X\t.\tstop_codon\t40\t42\t0\t-\t.\tgene_id D; transcript_id D1;\n" +
			"X\t.\tstart_codon\t91\t93\t0\t-\t.\tgene_id D; transcript_id D1;\n",
		GeneCnt:      1,
		FeatCnt:      1,
		IDs:          []string{"D"},
		Chrs:         []string{"X"},
		Starts:       []int{29},
		Ends:         []int{99},
		Orientations: []feat.Orientation{feat.Reverse},
	},
	{
		Name:  "Missing gene_id",
		Input: "X\t.\texon\t30\t50\t0\t-\t.\ttranscript_id D1;\n",
		Error: "gff: empty grouping gene_id field",
	},
	{
		Name:  "Missing transcript_id",
		Input: "X\t.\texon\t30\t50\t0\t-\t.\tgene_id D;\n",
		Error: "gff: empty grouping transcript_id field",
	},
	{
		Name: "Features on different location",
		Input: "" +
			"X\t.\texon\t30\t50\t0\t-\t.\tgene_id F; transcript_id F1;\n" +
			"Y\t.\texon\t80\t99\t0\t-\t.\tgene_id F; transcript_id F1;\n",
		Error: "geneio: features on varying location for gene F",
	},
	{
		Name: "Features on different orientation",
		Input: "" +
			"X\t.\texon\t30\t50\t0\t-\t.\tgene_id G; transcript_id G1;\n" +
			"X\t.\texon\t80\t99\t0\t+\t.\tgene_id G; transcript_id G1;\n",
		Error: "geneio: features with varying orientation for gene G",
	},
	{
		Name: "Single non coding gene with adjacent exons",
		Input: "" +
			"X\t.\texon\t30\t50\t0\t-\t.\tgene_id H; transcript_id H1;\n" +
			"X\t.\texon\t51\t99\t0\t-\t.\tgene_id H; transcript_id H1;\n",
		GeneCnt:      1,
		FeatCnt:      1,
		IDs:          []string{"H"},
		Chrs:         []string{"X"},
		Starts:       []int{29},
		Ends:         []int{99},
		Orientations: []feat.Orientation{feat.Reverse},
	},
	{
		Name: "Single non coding gene with overlapping exons",
		Input: "" +
			"X\t.\texon\t2\t70\t0\t+\t.\tgene_id F; transcript_id F1;\n" +
			"X\t.\texon\t70\t90\t0\t+\t.\tgene_id F; transcript_id F1;\n",
		Error: "exons overlap",
	},
	{
		Name: "Single coding gene with adjacent exons",
		Input: "" +
			"X\t.\texon\t10\t100\t0\t-\t.\tgene_id I; transcript_id I1;\n" +
			"X\t.\tstop_codon\t40\t42\t0\t-\t.\tgene_id I; transcript_id I1;\n" +
			"X\t.\tstart_codon\t91\t93\t0\t-\t.\tgene_id I; transcript_id I1;\n",
		GeneCnt:      1,
		FeatCnt:      1,
		IDs:          []string{"I"},
		Chrs:         []string{"X"},
		Starts:       []int{9},
		Ends:         []int{100},
		Orientations: []feat.Orientation{feat.Reverse},
	},
	{
		Name: "Single coding gene with overlapping exons",
		Input: "" +
			"X\t.\texon\t2\t70\t0\t+\t.\tgene_id F; transcript_id F1;\n" +
			"X\t.\texon\t70\t90\t0\t+\t.\tgene_id F; transcript_id F1;\n" +
			"X\t.\tstart_codon\t60\t62\t0\t+\t.\tgene_id F; transcript_id F1;\n" +
			"X\t.\tstop_codon\t81\t83\t0\t+\t.\tgene_id F; transcript_id F1;\n",
		Error: "exons overlap",
	},
	{
		Name: "Alternative group tags",
		Input: "" +
			"X\t.\texon\t10\t20\t0\t+\t.\tgid A; tid A1;\n" +
			"X\t.\texon\t15\t30\t0\t+\t.\tgid A; tid A2;\n" +
			"X\t.\texon\t20\t45\t0\t+\t.\tgid A; tid A3;\n" +
			"X\t.\texon\t25\t30\t0\t+\t.\tgid A; tid A4;\n" +
			"Y\t.\texon\t50\t90\t0\t-\t.\tgid B; tid B1;\n" +
			"Y\t.\tfoo\t100\t200\t0\t-\t.\tgid B; tid B1;\n" +
			"Y\t.\tstop_codon\t60\t62\t0\t-\t.\tgid B; tid B1;\n" +
			"Y\t.\tstart_codon\t71\t73\t0\t-\t.\tgid B; tid B1;\n",
		GeneTag:       "gid",
		TranscriptTag: "tid",
		GeneCnt:       2,
		FeatCnt:       5,
		IDs:           []string{"A", "B"},
		Chrs:          []string{"X", "Y"},
		Starts:        []int{9, 49},
		Ends:          []int{45, 90},
		Orientations:  []feat.Orientation{feat.Forward, feat.Reverse},
	},
}

func TestRead(t *testing.T) {
	for _, tt := range readTests {
		r := NewReader(gff.NewReader(strings.NewReader(tt.Input)))
		if tt.GeneTag != "" {
			_ = r.SetGeneTag(tt.GeneTag)
		}
		if tt.TranscriptTag != "" {
			_ = r.SetTranscriptTag(tt.TranscriptTag)
		}

		genes, err := r.ReadAll()
		if tt.Error != "" {
			if err == nil || !strings.Contains(err.Error(), tt.Error) {
				t.Errorf("%s: error %q, want error %q", tt.Name, err, tt.Error)
			}
			continue
		} else if err != nil {
			t.Errorf("%s: unexpected error %q", tt.Name, err)
			continue
		}

		if err := r.SetGeneTag("foo"); err == nil {
			t.Errorf("%s: expected error %q", tt.Name,
				"gff: cannot set GeneTag after first call to Read")
		}
		if err := r.SetTranscriptTag("foo"); err == nil {
			t.Errorf("%s: expected error %q", tt.Name,
				"gff: cannot set TranscriptTag after first call to Read")
		}

		var feats []feat.Feature
		var ids, chrs []string
		var starts, ends []int
		var orientations []feat.Orientation
		for _, g := range genes {
			feats = append(feats, g.Features()...)
			ids = append(ids, g.Name())
			starts = append(starts, g.Start())
			ends = append(ends, g.End())
			chrs = append(chrs, g.Location().Name())
			orientations = append(orientations, g.Orientation())
		}
		if len(genes) != tt.GeneCnt {
			t.Errorf("%s: out gene count=%d want %d", tt.Name, len(genes), tt.GeneCnt)
		}
		if len(feats) != tt.FeatCnt {
			t.Errorf("%s: out feat count=%d want %d", tt.Name, len(feats), tt.FeatCnt)
		}
		if !reflect.DeepEqual(ids, tt.IDs) {
			t.Errorf("%s: out ids=%q want %q", tt.Name, ids, tt.IDs)
		}
		if !reflect.DeepEqual(starts, tt.Starts) {
			t.Errorf("%s: out starts=%d want %d", tt.Name, starts, tt.Starts)
		}
		if !reflect.DeepEqual(ends, tt.Ends) {
			t.Errorf("%s: out ends=%d want %d", tt.Name, ends, tt.Ends)
		}
		if !reflect.DeepEqual(chrs, tt.Chrs) {
			t.Errorf("%s: out chrs=%q want %q", tt.Name, chrs, tt.Chrs)
		}
		if !reflect.DeepEqual(orientations, tt.Orientations) {
			t.Errorf("%s: out orientations=%q want %q", tt.Name, orientations, tt.Orientations)
		}
	}
}

func BenchmarkReadSorted(b *testing.B) {
	for i := 0; i < b.N; i++ {
		tt := readTests[0]
		r := NewReader(gff.NewReader(strings.NewReader(tt.Input)))
		_, _ = r.ReadAll()
	}
}

func Example() {
	data := `Y	.	exon	10	20	0	-	.	gene_id A; transcript_id A1
Y	.	exon	50	90	0	-	.	gene_id A; transcript_id A1;
Y	.	stop_codon	60	62	0	-	.	gene_id A; transcript_id A1
Y	.	start_codon	71	73	0	-	.	gene_id A; transcript_id A1`

	r := NewReader(gff.NewReader(strings.NewReader(data)))
	for {
		g, err := r.Read()
		if err != nil {
			if err == io.EOF {
				break
			} else {
				fmt.Println(err)
			}
		}
		fmt.Printf("gene: %s, transcript: %s", g.Name(), g.Features()[0].Name())

		// Output:
		// gene: A, transcript: A1
	}
}
