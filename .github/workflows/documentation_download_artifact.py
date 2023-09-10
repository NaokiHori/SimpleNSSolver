import sys
import datetime
import zipfile
import io
import requests

def set_params():
    # list all all artifacts (name and destination)
    artifacts = (
        {
            "name": "typical-cases-2d",
            "path": "docs/source/examples/typical/data/",
        },
        {
            "name": "typical-cases-3d",
            "path": "docs/source/examples/typical/data/",
        },
        {
            "name": "check-energy-conservations-2d",
            "path": "docs/source/examples/energy/data/",
        },
        {
            "name": "check-energy-conservations-3d",
            "path": "docs/source/examples/energy/data/",
        },
        {
            "name": "check-nusselt-agreements-1.e-1-2d",
            "path": "docs/source/examples/nu/data/",
        },
        {
            "name": "check-nusselt-agreements-1.e+0-2d",
            "path": "docs/source/examples/nu/data/",
        },
        {
            "name": "check-nusselt-agreements-1.e+1-2d",
            "path": "docs/source/examples/nu/data/",
        },
        {
            "name": "check-nusselt-agreements-1.e-1-3d",
            "path": "docs/source/examples/nu/data/",
        },
        {
            "name": "check-nusselt-agreements-1.e+0-3d",
            "path": "docs/source/examples/nu/data/",
        },
        {
            "name": "check-nusselt-agreements-1.e+1-3d",
            "path": "docs/source/examples/nu/data/",
        },
        {
            "name": "gl",
            "path": "docs/source/examples/gl/data/",
        },
    )
    return artifacts

def prepare_curl_header(token):
    headers = {
            "Accept": "application/vnd.github+json",
            "Authorization": f"Bearer {token}",
            "X-GitHub-Api-Version": "2022-11-28",
    }
    return headers

def kernel_get(url, headers):
    r = requests.get(
            url=url,
            headers=headers
    )
    if not r.ok:
        msg = f"request failed: {url}"
        raise RuntimeError(msg)
    return r

def list_artifacts(token):
    r = kernel_get(
            url="https://api.github.com/repos/NaokiHori/SimpleNSSolver/actions/artifacts",
            headers=prepare_curl_header(token)
    )
    artifacts = r.json()
    return artifacts

def get_and_save_artifact(token, param, artifacts):
    def convert_time(strtime):
        # from string to datetime.datetime
        # string: e.g. '2023-04-02T11:50:06Z'
        keywords = dict()
        # year
        year = strtime[:4]
        keywords["year"] = int(year)
        strtime = strtime.removeprefix(f"{year}-")
        # month
        month = strtime[:2]
        keywords["month"] = int(month)
        strtime = strtime.removeprefix(f"{month}-")
        # day
        day = strtime[:2]
        keywords["day"] = int(day)
        strtime = strtime.removeprefix(f"{day}T")
        # hour
        hour = strtime[:2]
        keywords["hour"] = int(hour)
        strtime = strtime.removeprefix(f"{hour}:")
        # minute
        minute = strtime[:2]
        keywords["minute"] = int(minute)
        strtime = strtime.removeprefix(f"{minute}:")
        # second
        second = strtime[:2]
        keywords["second"] = int(second)
        return datetime.datetime(**keywords)
    # initialise with N/A data
    result = {
            "url": "",
            "time": datetime.datetime(year=datetime.MINYEAR, month=1, day=1)
    }
    # check all artifacts and take the latest one
    for a in artifacts["artifacts"]:
        if a["name"] == param["name"]:
            time = a["created_at"]
            time = convert_time(time)
            if time > result["time"]:
                result["url"] = a["archive_download_url"]
                result["time"] = time
    if None == result["url"]:
        msg = "artifact '{}' is not found".format(param["name"])
        raise RuntimeError(msg)
    # get and save to the desired place
    r = kernel_get(
            url=result["url"],
            headers=prepare_curl_header(token)
    )
    zipfile.ZipFile(io.BytesIO(r.content)).extractall(param["path"])

def main(token):
    params = set_params()
    artifacts = list_artifacts(token)
    for param in params:
        get_and_save_artifact(token, param, artifacts)

if __name__ == "__main__":
    argv = sys.argv
    assert(2 == len(argv))
    main(argv[1])

