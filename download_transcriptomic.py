import requests
from bs4 import BeautifulSoup
from urllib.parse import urlparse, parse_qs

def fetch_download_links(base_url):
    try:
        response = requests.get(base_url)
        response.raise_for_status()
        soup = BeautifulSoup(response.text, 'lxml')
        links = soup.find_all('a', href=True)
        download_links = []

        for link in links:
            href = link['href']
            parsed = urlparse(href)
            # 检查路径是否为/downloadData/deDownload
            if parsed.path == "/downloadData/deDownload":
                query_params = parse_qs(parsed.query)
                path_value = query_params.get("path", [""])[0]
                # 检查path参数是否以.fa.gz结尾
                if path_value.endswith(".fa.gz"):
                    full_url = f"http://mgbase.qnlm.ac{href}"
                    download_links.append(full_url)

        return download_links

    except requests.exceptions.RequestException as e:
        print(f"请求失败: {e}")
        return []

# 示例使用
base_url = "http://mgbase.qnlm.ac/page/download/deNovoDownload"
download_links = fetch_download_links(base_url)

if download_links:
    print("找到的下载链接如下：")
    for link in download_links:
        print(link)
else:
    print("未找到符合条件的下载链接。")
