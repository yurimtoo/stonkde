#!/usr/bin/env python

import sys, os
import copy
from datetime import datetime
import time
import re
import requests
import traceback
import json
import subprocess
import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as ss

plt.ion()


# Global variable, since this looks goofy when included in function kwargs
DELIM = "\n================================================================================\n"
 
def downloadFromUrl(filename='comments.txt', object_type='comment', 
                    update=True, delay=1.0, delim=DELIM):
    """
    Downloads all relevant comments.

    Inputs
    ------
    filename: str.  Output filename.
    object_type: str. Type of content to download.  
                      Options: 'comment' or 'submission'.
                      Default: 'comment'
    update: bool. Determines whether to update an existing `filename` that was previously downloaded.
    delay: float. Delay time in seconds between each request.  
                  Default: 1.00
    delim: str. Delimiter between each `object_type` saved to `filename`.  
                Default: a bunch of equals signs.

    Outputs
    -------
    Saves all `object_type` content into `filename`.

    Revisions
    -----
    xx/xx/xxxx  u/geppetto123   Original implementation
    22/12/2021  u/yurimtoo      Added delay/delim variables and related code, 
                                as well as the update functionality
    23/12/2021  u/yurimtoo      Refactored code into a sequence of functions 
                                rather than scripts
    """
    # Sanity checks
    assert object_type in ['comment', 'submission']
    assert delay > 0

    # Set necessary parameters
    username = "Roid_Rage_Smurf"  # put the username you want to download in the quotes
    subreddit = ""  # put the subreddit you want to download in the quotes
    # leave either one blank to download an entire user's or subreddit's history
    # or fill in both to download a specific users history from a specific subreddit

    filter_string = None
    if username == "" and subreddit == "":
        print("Fill in either username or subreddit")
        sys.exit(0)
    elif username == "" and subreddit != "":
        filter_string = f"subreddit={subreddit}"
    elif username != "" and subreddit == "":
        filter_string = f"author={username}"
    else:
        filter_string = f"author={username}&subreddit={subreddit}"

    url = "https://api.pushshift.io/reddit/{}/search?limit=1000&sort=desc&{}&before="
    urlSubmission = "https://api.pushshift.io/reddit/search/submission/?ids={}&filter=parent_id,body,author,subreddit,created_utc,id"

    start_time = datetime.utcnow()

    # If updating an existing file, get the last content downloaded when script last executed
    if update:
        assert os.path.exists(filename)
        with open(filename, 'r') as foo:
            foostr = foo.read()
        # Get the last comment downloaded when script last executed
        newest = foostr.split(delim, 1)[0]
        # Adjust the saved filename, so that later they can be combined
        orig_foo = copy.copy(filename)
        filename = filename.replace('.', '_update.')
    print(f"Saving {object_type}s to {filename}")
 
    count = 0
    handle = open(filename, 'w')
    previous_epoch = int(start_time.timestamp())
    while True:
        new_url = url.format(object_type, filter_string)+str(previous_epoch)
        json_text = requests.get(new_url, headers={'User-Agent': "Post downloader by /u/Watchful1"})
        time.sleep(delay)  # pushshift has a rate limit, if we send requests too fast it will start returning error messages
        try:
            json_data = json_text.json()#
            #print (json_data)
        except json.decoder.JSONDecodeError:
            time.sleep(2)
            continue
 
        if 'data' not in json_data:
            break
        objects = json_data['data']
        if len(objects) == 0:
            break
 
        for obj in objects:
            previous_epoch = obj['created_utc'] - 1
            count += 1
            if object_type == 'comment':
                body = obj['body'].encode(encoding='ascii', errors='ignore').decode()
                if update:
                    if body in newest:
                        handle.close()
                        merge(orig_foo, filename)
                        return
                try:
                    handle.write(f"https://www.reddit.com{obj['permalink']}\n")
                    handle.write("This was created on: "+datetime.fromtimestamp(obj['created_utc']).strftime("%d.%m.%Y %H:%M:%S") + "\n")
                    # Get Thread username with max. 5 attemps
                    usernameattempts = 0
                    while usernameattempts < 5:
                        newurlsubmission = urlSubmission.format(obj['link_id'])
                        json_textsubmission = requests.get(newurlsubmission, headers={'User-Agent': "Post downloader by /u/Watchful1"})
                        time.sleep(delay * (usernameattempts+1))
                        handle.write(f"Getting username...\n")
                        try:
                            json_datasubmission = json_textsubmission.json()
                        except json.decoder.JSONDecodeError:
                            time.sleep(2)
                            usernameattempts += 1
                            continue
 
                        if 'data' not in json_data or len(json_datasubmission['data'])==0:
                            usernameattempts=6
 
                        if usernameattempts < 5:
                            userobject = json_datasubmission['data']
                            threadusername = userobject[0]['author']
                            handle.write(f"Thread by User: /u/{threadusername}\n")
                            print("Current Count {} with Thread User: {} in {}".format(count,threadusername,datetime.fromtimestamp(previous_epoch).strftime("%d.%m.%Y %H:%M:%S")))
                            break
                        else:
                            handle.write(f"Couldn't get parentuser for comment: https://www.reddit.com{obj['permalink']}\n")
 
                    # Print Body
                    handle.write(body)
                    handle.write(delim)
                    handle.flush()
                    os.fsync(handle)
                except Exception as err:
                    print(f"Couldn't print comment: https://www.reddit.com{obj['permalink']}")
                    print(traceback.format_exc())
            elif object_type == 'submission':
                if obj['is_self']:
                    if 'selftext' not in obj:
                        continue
                    try:
                        body = obj['selftext'].encode(encoding='ascii', errors='ignore').decode()
                        if update:
                            if body in newest:
                                handle.close()
                                merge(orig_foo, filename)
                                return
                        handle.write(str(obj['score']))
                        handle.write(" : ")
                        handle.write(datetime.fromtimestamp(obj['created_utc']).strftime("%Y-%m-%d"))
                        handle.write("\n")
                        handle.write(body)
                        handle.write(delim)
                        handle.flush()
                        os.fsync(handle)
                    except Exception as err:
                        print(f"Couldn't print post: {obj['url']}")
                        print(traceback.format_exc())
 
        print("Saved {} {}s through {}".format(count, object_type, datetime.fromtimestamp(previous_epoch).strftime("%Y-%m-%d")))
 
    print(f"Saved {count} {object_type}s")
    handle.close()


def merge(oldfoo, newfoo):
    """
    Merges the newly-downloaded data with previously-downloaded data.

    Inputs
    ------
    oldfoo: str. Previously-downloaded data.
    newfoo: str. Newly-downloaded data.

    Revisions
    ---------
    23/12/2021  u/yurimtoo      Original implementation.
    """
    subprocess.run('cat '+oldfoo+' >> '+newfoo, shell=True)
    subprocess.run(['mv', newfoo, oldfoo])
    return


def reduce_data(filename="comments.txt"):
    """
    Reduces the downloaded data to the number of shares held by each unique user

    Inputs
    ------
    filename: str. File that contains the downloaded data.

    Outputs
    -------
    unique_nshares: array, float.  Share counts held by each unique user 
                                   that has fed the DRS bot.

    Revisions
    ---------
    xx/xx/xxxx  u/geppetto123   Original implementation
    14/12/2021  u/yurimtoo      Implemented KDE, updated plotting
    23/12/2021  u/yurimtoo      Refactored to work with the raw downloaded data
    24/12/2012  u/yurimtoo      Added logic to remove duplicates
    """
    # parsing the comments.txt file
    textfile = open(filename, 'r')
    filetext = textfile.read()
    textfile.close()
     
    matches = re.findall("===============================\n(https:\/\/www\.reddit\.com\/.*?)\n[\s\S]{0,50}(\d\d.\d\d.\d\d\d\d) (\d\d:\d\d:\d\d)[\s\S]{0,100}(\/u\/.*)\n[\s\S]{5,300}[^\d](\d*) shares added", filetext, flags=re.IGNORECASE)
    resets = re.findall("===============================\n(https:\/\/www\.reddit\.com\/.*?)\n[\s\S]{0,50}(\d\d.\d\d.\d\d\d\d) (\d\d:\d\d:\d\d)[\s\S]{0,100}(\/u\/.*)\n[\s\S]{5,300}[^\d] ALL ENTRIES ARE RESET", filetext, flags=re.IGNORECASE)

    # Each match is link, date, time, username, shares
    # We only need the username and shares
    links   = []
    dates   = []
    times   = []
    users   = []
    nshares = []
    for singlematch in matches:
        if singlematch[3] == '[deleted]':
            continue
        try:
            links.append(singlematch[0])
            dates.append(singlematch[1])
            times.append(singlematch[2])
            users.append(singlematch[3])
            nshares.append(singlematch[4])
        except Exception as err:
            print(f"Couldn't print post:{singlematch[1]};{singlematch[2]};{singlematch[3]};{singlematch[4]};{singlematch[0]}")
            print(traceback.format_exc())

    # Sanity check
    assert len(users) == len(nshares)

    # Get set of unique users with each of their posts
    unique_users = []
    user_links = []
    user_dates = []
    user_times = []
    user_share = []
    for i in range(len(users)):
        # Some comments are from the user that runs the bot, don't store those
        if users[i] == '/u/Roid_Rage_Smurf':
            continue
        # Someone tricked the bot early on with an absurd sharecount, omit those
        if int(nshares[i]) > 1e6:
            continue
        try:
            ind = unique_users.index(users[i])
            user_links[ind].append(links[i])
            user_dates[ind].append(dates[i])
            user_times[ind].append(times[i])
            user_share[ind].append(nshares[i])
        except:
            unique_users.append(users[i])
            user_links.append([links[i]])
            user_dates.append([dates[i]])
            user_times.append([times[i]])
            user_share.append([nshares[i]])

    # Now reduce this data and remove duplicates
    for i in range(len(unique_users)):
        if len(user_links) == 1:
            # Only 1 post DRSing, so can't have duplicates
            continue
        these_shares = np.asarray(user_share[i])
        if np.unique(these_shares).size == these_shares.size:
            # Each post DRS'd different number of shares, so can't have duplicates
            continue
        else:
            badinds = []
            # Some posts have the same number of shares DRS'd. Verify that they are not duplicates
            for j in range(these_shares.size):
                for k in range(j+1, these_shares.size):
                    if these_shares[j] == these_shares[k]:
                        if 'GMEOrphans' in user_links[i][j]:
                            # All DRS'd shares are in the same thread, so check using the time
                            date1 = [int(val) for val in user_dates[i][j].split('.')]
                            date2 = [int(val) for val in user_dates[i][k].split('.')]
                            time1 = [int(val) for val in user_times[i][j].split(':')]
                            time2 = [int(val) for val in user_times[i][k].split(':')]
                            a = datetime(date1[2], date1[1], date1[0], time1[0], time1[1], time1[2])
                            b = datetime(date2[2], date2[1], date2[0], time2[0], time2[1], time2[2])
                            c = a - b
                            if c.total_seconds()/60. < 10:
                                # Comments within 10 minutes are almost surely duplicates
                                badinds.append(k)
                        else:
                            # For any sub besides r/GMEOrphans, it is 1 thread per set of shares DRS'd
                            if user_links[i][j].rsplit('/', 2)[0] == user_links[i][k].rsplit('/', 2)[0]:
                                # Same thread -- must be duplicate (e.g., comment deleted and re-posted)
                                badinds.append(k)
            # Remove the duplicates
            badinds = sorted(np.unique(np.asarray(badinds)), reverse=True)
            for idx in badinds:
                user_links[i].pop(idx)
                user_dates[i].pop(idx)
                user_times[i].pop(idx)
                user_share[i].pop(idx)

    # Process the DRS reset posts
    reset_link = []
    reset_date = []
    reset_time = []
    reset_user = []
    for reset in resets:
        reset_link.append(reset[0])
        reset_date.append(reset[1])
        reset_time.append(reset[2])
        reset_user.append(reset[3])

    # Remove all data after a user submitted a reset request
    for i in range(len(reset_user)):
        try:
            iusr = unique_users.index(reset_user[i])
            # Build datetime object for the reset request
            thisday   = [int(val) for val in reset_date[i].split('.')]
            thistime  = [int(val) for val in reset_time[i].split(':')]
            thisreset = datetime(thisday [2], thisday [1], thisday [0], 
                                 thistime[0], thistime[1], thistime[2])
            # Compare this datetime to each DRS post, 
            # starting from the first post (work backwards thru list)
            theseinds = np.linspace(len(user_dates[iusr])-1, 0, len(user_dates[iusr]), dtype=int)
            for j in theseinds:
                oldday  = [int(val) for val in user_dates[iusr][j].split('.')]
                oldtime = [int(val) for val in user_times[iusr][j].split(':')]
                if thisreset > datetime(oldday [2], oldday [1], oldday [0], 
                                        oldtime[0], oldtime[1], oldtime[2]):
                    # Reset is more recent, remove this instance
                    user_links[iusr].pop(j)
                    user_dates[iusr].pop(j)
                    user_times[iusr].pop(j)
                    user_share[iusr].pop(j)
                    if j==0 and len(user_dates[iusr]) == 0:
                        # Looks like a user hasn't fed the bot since they reset...
                        # Remove them
                        user_links.pop(iusr)
                        user_dates.pop(iusr)
                        user_times.pop(iusr)
                        user_share.pop(iusr)
                        unique_users.pop(iusr)
        except:
            # User does not have any DRS posts, only the reset request
            continue

    # Combine separate posts into total # of DRS'd shares for each user
    unique_nshares = np.zeros(len(unique_users))
    for i in range(len(unique_users)):
        unique_nshares[i] = np.sum(np.asarray(user_share[i], dtype=int))

    return unique_nshares


def plot(nshares):
    """
    Produces plots of the sharecounts.

    Inputs
    ------
    fout: str. Basename for output plots.
    nshares: array, float.  Share counts held by each unique user 
                            that has fed the DRS bot.
    nbins: int. Number of bins to use when plotting histogram.
    logbins: bool.  Determines whether to use logarithmically-spaced bins (True) 
                    or linearly-spaced bins (False).
    kde_bw: str or scalar.  Method to determine KDE bandwidth.
                  Options: 'scott', 'silvermann', or scalar constant
                  Default: 'scott'

    Outputs
    -------
    Various plots with the `fout` basename.

    Revisions
    ---------
    14/12/2021  u/yurimtoo      Original implementation
    23/12/2021  u/yurimtoo      Refactored into a function, 
                                added calls to archive() and make_gif()
    """
    # Plot the histogram with various bin sizes and spacings
    #bins = np.logspace(0, np.log10(nshares.max()), nbins)
    for nbins in [int(np.ceil(nshares.max())),20000,10000,5000,1000,100]:
        bins = np.linspace(1, int(np.ceil(nshares.max())), nbins)
        ncount, bins, patches = plt.hist(nshares, bins=bins, 
                                         density=True, 
                                         histtype='step', 
                                         label='Histogram, '+str(nbins)+' bins')
    plt.xlabel('Number of shares')
    #plt.xscale('log')
    #plt.yscale('log')
    plt.ylabel('Probability density')
    plt.legend(loc='best')
    plt.title('Distribution of $GME Shareholders\n(n='+str(nshares.size)+')')
    plt.gca().set_rasterized(True)
    plt.savefig('post/figs/histogram_binsize-comp_full.png', dpi=2000)
    xlims = plt.xlim()
    plt.xlim(-10, 1010)
    plt.savefig('post/figs/histogram_binsize-comp_under1k.png', dpi=2000)
    plt.yscale('log')
    plt.savefig('post/figs/histogram_binsize-comp_under1k_ylog.png', dpi=2000)
    plt.xlim(*xlims)
    plt.savefig('post/figs/histogram_binsize-comp_full_ylog.png', dpi=2000)
    plt.close()

    # 10,000 bins looks to do the trick
    nbins = 10000
    bins = np.linspace(1, nshares.max(), nbins)
    ncount, bins, patches = plt.hist(nshares, bins=bins, 
                                     density=True, 
                                     histtype='step', 
                                     label='Histogram, '+str(nbins)+' bins', zorder=40)

    # KDE to determine the underlying dist'n
    # Bandwidth selection is arguably the most important part of KDE
    # Scott's rule and Silverman's rule are generally good choices, 
    # particularly for higher dimensional data
    # But the DRS bot data is sparse, so those rules tend to over-smooth the data
    # The definition below matches Scott's rule when d=4, 
    # and it doesn't over-smooth as much as the aforementioned rules
    ape_bw      = nshares.size**(-1. / (len(nshares.shape) * 2.0))
    ape_bw_alt1 = nshares.size**(-1. / (len(nshares.shape) * 1.5))
    ape_bw_alt2 = nshares.size**(-1. / (len(nshares.shape) * 3.0))
    # Now do KDE for each method
    #rep = ss.gaussian_kde(nshares, bw_method=kde_bw)
    kde_scott    = ss.gaussian_kde(nshares, bw_method='scott')
    kde_silver   = ss.gaussian_kde(nshares, bw_method='silverman')
    kde_ape      = ss.gaussian_kde(nshares, bw_method=ape_bw)
    kde_ape_alt1 = ss.gaussian_kde(nshares, bw_method=ape_bw_alt1)
    kde_ape_alt2 = ss.gaussian_kde(nshares, bw_method=ape_bw_alt2)
    # Evaluate the KDE on a fine grid, and ensure it integrates to 1
    finegrid = np.linspace(1, nshares.max(), 100*nbins)
    #val = rep.evaluate(finegrid) # Evaluate the KDE according to the grid
    #val /= np.trapz(val, finegrid)
    scott  = kde_scott.evaluate(finegrid)
    scott /= np.trapz(scott, finegrid)
    silver  = kde_silver.evaluate(finegrid)
    silver /= np.trapz(silver, finegrid)
    ape  = kde_ape.evaluate(finegrid)
    ape /= np.trapz(ape, finegrid)
    ape_alt1 = kde_ape_alt1.evaluate(finegrid)
    ape_alt1 /= np.trapz(ape_alt1, finegrid)
    ape_alt2 = kde_ape_alt2.evaluate(finegrid)
    ape_alt2 /= np.trapz(ape_alt2, finegrid)

    plt.plot(finegrid, scott, label="KDE, Scott's Rule", zorder=40)
    plt.plot(finegrid, silver, label="KDE, Silverman's Rule", zorder=40)
    plt.plot(finegrid, ape, label="KDE, modified Scott's Rule", zorder=40)
    plt.plot(finegrid, ape_alt1, label="KDE, modified Scott's Rule, alt 1", ls='--', zorder=1)
    plt.plot(finegrid, ape_alt2, label="KDE, modified Scott's Rule, alt 2", ls='--', zorder=10)
    plt.legend(loc='best')
    plt.gca().set_rasterized(True)
    plt.savefig('post/figs/histogram_with_kde.png', dpi=1500)
    # Store current xlims
    xlims = plt.xlim()
    # Trim to under 1000 shares for prettier plots
    plt.xlim(-10, 1010)
    plt.savefig('post/figs/histogram_with_kde_under1k.png', dpi=1500)
    plt.ylim(1e-5, 1e-1)
    plt.xscale('log')
    plt.yscale('log')
    plt.xlim(1, xlims[1]) # back to orignial scale
    plt.savefig('post/figs/histogram_with_kde_log.png', dpi=1500)
    plt.close()
    # Let's look at the whales
    ncount, bins, patches = plt.hist(nshares, bins=bins, 
                                     density=True, 
                                     histtype='step', 
                                     label='Histogram, '+str(nbins)+' bins', zorder=1)
    plt.plot(finegrid, scott, label="KDE, Scott's Rule", zorder=100)
    plt.plot(finegrid, silver, label="KDE, Silverman's Rule", zorder=100)
    plt.plot(finegrid, ape, label="KDE, modified Scott's Rule", zorder=1)
    plt.xscale('log')
    plt.yscale('log')
    plt.xlim(390, 40100)
    plt.ylim(1e-8, 2e-3)
    plt.legend(loc='best')
    plt.gca().set_rasterized(True)
    plt.savefig('post/figs/histogram_with_kde_log_400-40k.png', dpi=1000)
    plt.close()
    # Archive these plots for the future
    archive()
    # Update the gif of the data
    make_gifs()

    return


def archive():
    """
    Archives the latest plots

    Revisions
    ---------
    23/12/2021  u/yurimtoo      Original implementation
    """
    pass


def make_gifs():
    """
    Create gifs of the plots

    Revisions
    ---------
    23/12/2021  u/yurimtoo      Original implementation
    """
    pass


def funstats(nshares):
    print("Number of apes that have fed the DRS bot:", nshares.size)
    print("The median shares held is", np.median(nshares))
    print("The biggest whale has", nshares.max(), "shares, which is " + str(100*nshares.max()/nshares.sum())[:5] + "% of all DRS'd shares")
    iX    = nshares<10
    iXX   = (nshares<100)  * (nshares >=10)
    iXXX  = (nshares<1000) * (nshares >=100)
    iXXXX = nshares>=1000
    print(str(nshares[nshares==1].size/nshares.size*100)[:6]+"% of holders have 1 share")
    print(str(nshares[iX   ].size/nshares.size*100)[:6]+"% of holders are X holders")
    print(str(nshares[iXX  ].size/nshares.size*100)[:6]+"% of holders are XX holders")
    print(str(nshares[iXXX ].size/nshares.size*100)[:6]+"% of holders are XXX holders")
    print(str(nshares[iXXXX].size/nshares.size*100)[:6]+"% of holders are XXXX+ holders")
    print(str(nshares[nshares<=   2].size/nshares.size*100)[:6]+"% of holders have <=  2 shares")
    print(str(nshares[nshares<=   4].size/nshares.size*100)[:6]+"% of holders have <=  4 shares")
    print(str(nshares[nshares<=   8].size/nshares.size*100)[:6]+"% of holders have <=  8 shares")
    print(str(nshares[nshares<=  16].size/nshares.size*100)[:6]+"% of holders have <= 16 shares")
    print(str(nshares[nshares<=  32].size/nshares.size*100)[:6]+"% of holders have <= 32 shares")
    print(str(nshares[nshares<=  64].size/nshares.size*100)[:6]+"% of holders have <= 64 shares")
    print(str(nshares[nshares<= 128].size/nshares.size*100)[:6]+"% of holders have <=128 shares")
    # What portion of total DRS'd shares?
    print("X     apes hold " + str(100*(nshares[iX].sum())/nshares.sum())[:5] + "% of DRS'd shares")
    print("XX    apes hold " + str(100*(nshares[iXX].sum())/nshares.sum())[:5] + "% of DRS'd shares")
    print("XXX   apes hold " + str(100*(nshares[iXXX].sum())/nshares.sum())[:5] + "% of DRS'd shares")
    print("XXXX+ apes hold " + str(100*(nshares[iXXXX].sum())/nshares.sum())[:5] + "% of DRS'd shares")
    # Find the sharecount that splits the total DRS'd shares into two nearly equal groups
    sortedshares = np.sort(nshares) #not to be confused with shorted shares...
    cumsum = np.cumsum(sortedshares)
    ncut = sortedshares[np.argmin(np.abs(cumsum - cumsum[-1]/2.))]
    print("Apes holding", ncut, "or fewer shares hold a total of ~50% of DRS'd shares")
    print("Quantiles:")
    print("  10%:", np.quantile(nshares, 0.10), "shares")
    print("  25%:", np.quantile(nshares, 0.25), "shares")
    print("  33%:", np.quantile(nshares, 0.33), "shares")
    print("  69%:", np.quantile(nshares, 0.69), "shares")
    print("  75%:", np.quantile(nshares, 0.75), "shares")
    print("  90%:", np.quantile(nshares, 0.90), "shares")
    print("  95%:", np.quantile(nshares, 0.95), "shares")
    print("  99%:", np.quantile(nshares, 0.99), "shares")
    print("Sharecounts of the 20 biggest whales that have fed the bot:")
    print(np.sort(nshares)[::-1][:20])


if __name__ == '__main__':
    # Download comments
    downloadFromUrl("comments.txt", "comment", update=True)
    # Reduce the data and plot
    nshares = reduce_data("comments.txt")
    funstats(nshares)
    plot(nshares)


