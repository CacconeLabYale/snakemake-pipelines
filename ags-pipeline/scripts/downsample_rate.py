# import ipdb

import bokeh as bk
from bokeh import mpl
from bokeh.plotting import figure as f
import bokeh.charts as c

import os
import glob
import pandas as pd



def get_sample_rate(original, sample_to):
    if original <= sample_to:
        return 1

    else:
        return sample_to/original

def add_dwnsmpl_rate(df):
    sample_to = df.mapped_reads.min() * 2

    df['sample_rate'] = df.mapped_reads.apply(lambda x: get_sample_rate(x, sample_to))
    df['sampled_read_coverage'] = df['sample_rate'] * df.mapped_reads



def main(count_file, fig_dir, read_count_files=None, file_pattern=None):

    headers = ["individual","mapped_reads"]

    if read_count_files is None:
        read_count_files = glob.glob('{file_pattern}'.format(file_pattern=file_pattern))

    read_count_data = []

    for path in read_count_files:
            bam = os.path.splitext(os.path.basename(path))[0]
            count = open(path, 'rU').readline().rstrip('\n')
            read_count_data.append((bam, int(count)))


    # ipdb.set_trace()
    mapped_reads_df = pd.DataFrame(data=read_count_data, index=None, columns=headers)
    mapped_reads_df['all'] = 'all'

    # # plot read numbers and distribution
    # box = c.BoxPlot(data=mapped_reads_df, label="all", values="mapped_reads",
    #               color=None, group=None,
    #               xscale='categorical', yscale='linear', xgrid=False, ygrid=True, continuous_range=None)
    # 
    # 
    # bar = c.Bar(data=mapped_reads_df, label="file", values="mapped_reads",
    #           color=None, stack=None, group=None,
    #           agg='sum',
    #           xscale='categorical', yscale='linear', xgrid=False, ygrid=True, continuous_range=None)
    # p = c.hplot(bar,box)
    # 
    # c.output_file('{fig_dir}/original_read_numbers.html'.format(fig_dir=fig_dir))
    # c.save(c.vplot(p))


    # # Calculate the percentage to downsample for each individual
    #
    # Based on disscussion from meeting on 2015-09-23, I will downsample all to within 2 fold excess of the lowest fly (MS11_0017 = 1 396 521 reads).

    add_dwnsmpl_rate(mapped_reads_df)

    # bar_old = c.Bar(data=mapped_reads_df, label="file", values="mapped_reads",
    #                 color=None, stack=None, group=None,
    #                 agg='sum',
    #                 xscale='categorical', yscale='linear', xgrid=False, ygrid=True, continuous_range=None)
    # 
    # 
    # bar_new  = c.Bar(data=mapped_reads_df, label="file", values="sampled_read_coverage",
    #                  color=None, stack=None, group=None,
    #                  agg='sum',
    #                  xscale='categorical', yscale='linear', xgrid=False, ygrid=True,
    #                  continuous_range=None)
    # 
    # bar_new.y_range = bar_old.y_range
    # 
    # c.output_file('{fig_dir}/compare_orig_to_dwnsampled.html'.format(fig_dir=fig_dir))
    # c.save(c.vplot(c.hplot(bar_old,bar_new)))

    # write out read number table with calculated downsample rates.
    mapped_reads_df.to_csv(path_or_buf=count_file, sep=',', index=False)

if __name__ == "__main__":
    main()
