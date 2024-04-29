import axios, { AxiosInstance } from 'axios';
/**
 * Interface for parameters required by GoogleCustomSearch class.
 */
export interface GoogleCustomSearchParams {
  apiKey?: string;
  googleCSEId?: string;
}

/**
 * Class that uses the Google Search API to perform custom searches.
 * Requires environment variables `GOOGLE_API_KEY` and `GOOGLE_CSE_ID` to
 * be set.
 */
export class GoogleCustomSearch {
  name = 'google-custom-search';

  protected apiKey: string;

  protected googleCSEId: string;

  protected instance: AxiosInstance;

  description =
    'a custom search engine. useful for when you need to answer questions about current events. input should be a search query. outputs a JSON array of results.';

  constructor(
    fields: GoogleCustomSearchParams = {
      apiKey: process.env.GOOGLE_SEARCH_API_KEY,
      googleCSEId: process.env.GOOGLE_SEARCH_API_CX,
    },
  ) {
    if (!fields.apiKey) {
      throw new Error(
        `Google API key not set. You can set it as "GOOGLE_API_KEY" in your environment variables.`,
      );
    }
    if (!fields.googleCSEId) {
      throw new Error(
        `Google custom search engine id not set. You can set it as "GOOGLE_CSE_ID" in your environment variables.`,
      );
    }
    this.apiKey = fields.apiKey;
    this.googleCSEId = fields.googleCSEId;
    const proxy =
      process.env.PROXY_ENABLE === 'true'
        ? {
            host: process.env.PROXY_HOST,
            port: +process.env.PROXY_PORT,
            protocol: process.env.PROXY_PROTOCOL,
          }
        : false;
    this.instance = axios.create({
      baseURL: 'https://www.googleapis.com/customsearch/v1',
      proxy,
    });
  }

  async call(input: string) {
    const response = await this.instance.get(
      `?key=${this.apiKey}&cx=${this.googleCSEId}&q=${encodeURIComponent(
        input,
      )}`,
    );
    if (response.data) {
      const json = response.data;
      const results =
        json?.items?.map(
          (item: { title?: string; link?: string; snippet?: string }) => ({
            title: item.title,
            link: item.link,
            snippet: item.snippet,
          }),
        ) ?? [];
      return JSON.stringify(results);
    } else {
      throw new Error(
        `Got ${response.status} error from Google custom search: ${response.statusText}`,
      );
    }
  }
}
