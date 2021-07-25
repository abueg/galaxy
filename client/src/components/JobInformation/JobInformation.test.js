import Vuex from "vuex";
import axios from "axios";
import MockAdapter from "axios-mock-adapter";
import { mount } from "@vue/test-utils";
import { getLocalVue } from "jest/helpers";
import JobInformation from "./JobInformation";
import datasetResponse from "components/DatasetInformation/testData/datasetResponse";
import jobResponse from "./testData/jobInformationResponse.json";

import flushPromises from "flush-promises";
import createCache from "vuex-cache";

jest.mock("app");
jest.mock("../History/caching");

const JOB_ID = "test_id";

const localVue = getLocalVue();

const jobStore = new Vuex.Store({
    plugins: [createCache()],
    modules: {
        jobStore: {
            actions: {
                fetchDataset: jest.fn(),
                fetchJob: jest.fn(),
            },
            getters: {
                dataset: (state) => (hda_id) => {
                    return datasetResponse;
                },
                job: (state) => (hda_id) => {
                    return jobResponse;
                },
            },
        },
    },
});

describe("JobInformation/JobInformation.vue", () => {
    let wrapper;
    let jobInfoTable;
    let axiosMock;

    beforeEach(() => {
        axiosMock = new MockAdapter(axios);
        axiosMock.onGet(new RegExp(`api/configuration/decode/*`)).reply(200, { decoded_id: 123 });
    });

    afterEach(() => {
        axiosMock.restore();
    });

    const verifyValues = (rendered_entries, infoTable, backendResponse) => {
        rendered_entries.forEach((entry) => {
            const renderedText = infoTable.find(`#${entry.id}`).text();
            const backendText = entry.wrapped_in_brackets
                ? `(${backendResponse[entry.backend_key]})`
                : backendResponse[entry.backend_key];
            // use includes because we append the decoded id
            expect(renderedText.includes(backendText.toString())).toBeTruthy();
        });
    };

    beforeEach(async () => {
        const propsData = {
            job_id: JOB_ID,
        };

        wrapper = mount(JobInformation, {
            store: jobStore,
            propsData,
            localVue,
        });
        jobInfoTable = wrapper.find("#job-information");
        await flushPromises();
    });

    it("job information table should be rendered", async () => {
        // table should exist
        expect(jobInfoTable).toBeTruthy();
        const rows = jobInfoTable.findAll("tr");
        // should contain 9 rows
        expect(rows.length).toBe(9);
    });

    it("stdout and stderr should be rendered", async () => {
        ["stdout", "stderr"].forEach((std) => {
            const label = jobInfoTable.find("#" + std);
            const value = label.find(".code");
            expect(value.text()).toBe(std);
        });
    });

    it("job messages", async () => {
        const rendered_link = jobInfoTable.findAll(`#job-messages li`);
        expect(rendered_link.length).toBe(jobResponse.job_messages.length);
        for (let i = 0; i < rendered_link.length; i++) {
            const msg = rendered_link.at(i).text();
            expect(jobResponse.job_messages.includes(msg));
        }
    });

    it("job_information API content", async () => {
        const rendered_entries = [
            { id: "galaxy-tool-id", backend_key: "tool_id" },
            { id: "galaxy-tool-version", backend_key: "tool_version" },
            { id: "encoded-job-id", backend_key: "id" },
            { id: "encoded-copied-from-job-id", backend_key: "copied_from_job_id" },
        ];
        verifyValues(rendered_entries, jobInfoTable, jobResponse);
    });
});
